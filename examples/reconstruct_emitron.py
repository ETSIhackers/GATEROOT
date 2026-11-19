"""Sequential GATEROOT prompts+delayeds -> preliminary OSEM -> SSS -> corrected OSEM.

Run with the existing EMITRON development environment, separately from the converter environment.
No activity truth, external normalization file, or scatter truth is supplied to reconstruction.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

CASES = [('IEC_6min',360,64,'umap_mIEC_ETSIPETscanner.hv'),
         ('Hoffman_3min',180,128,'umap_ETSIPETscanner256mmAFOV_brainSimplePhantom.hv'),
         ('Hoffman_6min',360,128,'umap_ETSIPETscanner256mmAFOV_brainSimplePhantom.hv')]


def save(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    return path


def make_config(source, out, duration, nz, attenuation, edges, resolution, cache, randoms=None):
    corrections={'attenuation':{'image':str(attenuation),'units':'cm^-1','projector':'siddon'}}
    if randoms:
        corrections.update(randoms=str(randoms),non_tof_distribution='uniform_time',tof_edges_ps=edges)
    return {'schema_version':'0.2',
        'image':{'n_voxels':[360,360,nz],'voxel_size_mm':[2.,2.,2.],'center_mm':[0.,0.,0.]},
        'frames_s':[[0.,float(duration)]],
        'acquisition':{'durations_s':[float(duration)],'starts_s':[0.],
                       'tracer':{'name':'unspecified','half_life_s':None},'decay_correction':'none'},
        'algorithm':'osem','iterations':2,'subsets':10,'backend':'cuda','projector':'siddon',
        'batch_size':100000,'output_prefix':str(out/'image'),'output_format':'both',
        'save':{'iteration_interval':0,'subset_interval':0,'always_final':True},
        'execution':{'device':'cuda:0','gpu_devices':[0],'cpu_threads':4,'backprojection_memory_mb':512},
        'iteration_cache':{'enabled':True,'directory':str(cache),'max_gib':12,
                           'decoded_events':True,'decoded_cache_directory':str(cache/'events')},
        'petsird_inspection':{'delayed_randoms_scale':1.,'delayed_randoms_mode':'histogram'},
        'tof':{'enabled':True,'resolution_fwhm_ps':resolution,
               'bin_width_ps':float(edges[1]-edges[0]),'continuous_bin_centers':True},
        'beds':[{'petsird_files':[str(source)],'offset_mm':[0.,0.,0.],
                 'exclude_lors_outside_fov':True,'corrections':corrections,
                 'sensitivity':{'source':'generate','tof_edges_ps':edges,
                    'generation':{'method':'exact','collapse_tof_bins':True,
                                  'checkpoint_file':str(out/'sensitivity.checkpoint.npz')}}}]}


def scatter_spec(out, edges):
    return {'method':'sss','activity_estimation':{'iterations':2},
        'output_file':str(out/'scatter.hscatter'),'batch_size':256,'max_memory_gib':2,
        'sampling':{'detectors_per_ring':16,'rings':8,'image_downsample_xyz':[5,5,1],
                    'lookup_batch_size':10000,'fit_batch_size':10000,
                    'max_sampled_pairs':100000,'max_fit_pairs':1000000000},
        'sss':{'backend':'cuda','device':'cuda:0','gpu_evaluation_batch':256,
               'max_ray_samples':1048576,'max_gpu_memory_gib':2,'progress':True,
               'ray_step_mm':4.,'scatter_stride':3,'max_pair_scatter_evaluations':500000000,
               'energy_window_keV':[430.,650.],'energy_resolution_fwhm_fraction':.11,
               'tail_margin_mm':8.,'tof_edges_ps':edges},
        'scaling':{'method':'tail_fit','min_bins':8}}


def run_case(source_dir, out_root, mu_dir, cache, case, units, plan):
    import numpy as np
    from emitron.petsird import read_header
    from emitron.projector import C_MM_PER_PS
    from emitron.image_io import read_image, write_image
    from emitron.config import ImageGrid
    name,duration,nz,mu_name=case
    source=(source_dir/(name+'_prompts_delayeds_norm.petsird')).resolve()
    mu_source=(mu_dir/mu_name).resolve()
    for path in (source,mu_source):
        if not path.is_file():raise FileNotFoundError(path)
    out=out_root/name
    out.mkdir(parents=True,exist_ok=True)
    header=read_header(source)
    scanner=header.scanner
    edges=(np.asarray(scanner.tof_bin_edges[0][0].edges,dtype=float)*2/C_MM_PER_PS).tolist()
    resolution=float(scanner.tof_resolution[0][0])*2/C_MM_PER_PS
    if not np.isfinite(resolution) or resolution<=0:raise ValueError('Invalid embedded TOF resolution')
    # These supplied phantom umaps were previously registered by centering their voxel grids
    # at the scanner origin. Preserve all values and spacing; do not flip or interpolate.
    attenuation=out/'mu511_centered.hdr'
    mu,meta=read_image(mu_source)
    if mu.ndim!=3 or not np.isfinite(mu).all() or np.any(mu<0):raise ValueError('Invalid 3-D umap')
    mu=np.asarray(mu,dtype=np.float32)*(10. if units=='mm^-1' else 1.)
    if any(n%f for n,f in zip(mu.shape[::-1],[5,5,1])):
        raise ValueError('This example requires x/y umap dimensions divisible by 5')
    grid=ImageGrid(tuple(mu.shape[::-1]),tuple(meta['voxel_size_mm']),(0.,0.,0.))
    signature={'source':str(source),'size':source.stat().st_size,'mtime_ns':source.stat().st_mtime_ns,
               'umap':str(mu_source),'umap_sha256':hashlib.sha256(mu.tobytes()).hexdigest(),
               'umap_units':units,'umap_center_mm':[0.,0.,0.],'duration_s':duration,'subsets':10}
    manifest=out/'example_inputs.json'
    if manifest.exists() and json.loads(manifest.read_text())!=signature:
        raise ValueError('Inputs changed: use a new output directory')
    if not attenuation.exists():write_image(attenuation.with_suffix(''),mu,grid,output_format='interfile')
    del mu
    if not manifest.exists():save(manifest,signature)
    cfg=make_config(source,out,duration,nz,attenuation,edges,resolution,cache)
    if plan:
        save(out/'plan.json',dict(configuration=cfg,scatter_estimation=scatter_spec(out,edges)))
        print(f'PLAN {name}: {source}; {duration}s, 2 iterations / 10 subsets; centered umap {attenuation}',flush=True)
        return
    if (out/'image_iter002.npy').exists():
        print(f'Retaining completed {name}',flush=True)
        return
    from emitron.settings import load_settings
    from emitron.recon import prepare_beds, reconstruct
    from emitron.execution import execution_scope
    from emitron.event_corrections import release_models
    # Generate a non-TOF histogram of measured delayed counts. Uniform timing is explicit
    # in the TOF model; the detector-pair distribution is not replaced by a global constant.
    descriptor=out/'randoms_descriptor.json'
    if descriptor.exists():
        randoms=Path(json.loads(descriptor.read_text())['randoms'])
        if not randoms.exists():raise FileNotFoundError(randoms)
    else:
        non_tof=deepcopy(cfg)
        non_tof['tof']={'enabled':False}
        non_tof['output_prefix']=str(out/'randoms')
        non_tof['beds'][0]['sensitivity']={'source':'generate','generation':{'method':'exact'}}
        path=save(out/'extract_delayeds.json',non_tof)
        settings,image_grid,frames=load_settings(path)
        with execution_scope(settings.get('execution')),release_models(settings):
            bed,_,_=prepare_beds(settings,image_grid,frames,with_sensitivity=False,
                                  audit_events=False,verify_prompts=False)[0]
            randoms=Path(bed['corrections']['randoms'])
        save(descriptor,{'randoms':str(randoms)})
    cfg=make_config(source,out,duration,nz,attenuation,edges,resolution,cache,randoms)
    scatter=out/'scatter.hscatter'
    sensitivity=out/'scatter_preliminary_bed000_sensitivity.npy'
    if not scatter.exists():
        from emitron.sss_interface import estimate_config
        cfg['beds'][0]['scatter_estimation']=scatter_spec(out,edges)
        # A restart after scatter fitting fails can reuse the already exported preliminary
        # sensitivity; it repeats only the two preliminary iterations and unfinished SSS.
        if sensitivity.exists():
            cfg['beds'][0].update(sensitivity={'source':'file','image':str(sensitivity),'tof_edges_ps':edges},
                sensitivity_includes_corrections=True,sensitivity_includes_petsird_efficiencies=True)
        estimate_config(save(out/'estimate_scatter.json',cfg))
    if not sensitivity.exists():raise FileNotFoundError(sensitivity)
    cfg['beds'][0].pop('scatter_estimation',None)
    cfg['beds'][0]['corrections']['scatter']=str(scatter)
    cfg['beds'][0].update(sensitivity={'source':'file','image':str(sensitivity),'tof_edges_ps':edges},
        sensitivity_includes_corrections=True,sensitivity_includes_petsird_efficiencies=True)
    final=out/'reconstruct.json'
    if not final.exists():
        from emitron.precompute_corrections import compile_config
        compile_config(save(out/'compile_additive.json',cfg),out/'additive.edh',final,
            storage='sidecar',interface='precorrected_additive',preserve_output_prefix=True)
    reconstruct(final)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mu-dir',type=Path,required=True)
    parser.add_argument('--mu-units',choices=['cm^-1','mm^-1'],required=True)
    parser.add_argument('--cache-dir',type=Path,required=True)
    parser.add_argument('--cases',nargs='+',choices=[c[0] for c in CASES],default=[c[0] for c in CASES])
    parser.add_argument('--plan',action='store_true',help='Check inputs/header/umap and write plans; no event scan or reconstruction')
    args=parser.parse_args()
    # Fail immediately on unsupported EMITRON versions, before processing huge event files.
    from inspect import signature
    from emitron.recon import prepare_beds
    if 'audit_events' not in signature(prepare_beds).parameters:
        raise RuntimeError('Use the EMITRON development checkout with validated PETSIRD event caching and SSS')
    args.cache_dir=args.cache_dir.resolve()
    args.cache_dir.mkdir(parents=True,exist_ok=True)
    for case in CASES:
        if case[0] in args.cases:
            run_case(args.input_dir.resolve(),args.output.resolve(),args.mu_dir.resolve(),
                     args.cache_dir,case,args.mu_units,args.plan)


if __name__=='__main__':main()
