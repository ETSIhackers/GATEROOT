"""Convert cylindrical ETSIPET GATE ROOT data to binary PETSIRD 0.13."""

import argparse
import json
from itertools import product
from pathlib import Path

from importlib.metadata import version
import numpy as np
import petsird as p
from petsird.helpers.geometry import construct_cuboid_box
import uproot


C_MM_PER_PS = .299792458
NCF = np.dtype([('factor','<f4'),('det1','<u4'),('det2','<u4')])
LEVELS = ('mod','smod','cry')
# coincidence categories from the GATE truth: same decay (eventID), and number of Compton interactions
# in the phantom for the two photons together (0: true, 1: single scatter, 2 or more: multiple scatter)
CATEGORIES = ('true','single-scatter','multiple-scatter','random')
SELECTION_GROUPS = {'all':CATEGORIES,'scatter':('single-scatter','multiple-scatter')}


def selection(text):
    """comma-separated coincidence categories (true, single-scatter, multiple-scatter, random) or groups (scatter, all)"""
    chosen=set()
    for word in text.split(','):
        if word not in CATEGORIES and word not in SELECTION_GROUPS:
            raise argparse.ArgumentTypeError(f'unknown selection {word!r}')
        chosen.update(SELECTION_GROUPS.get(word,(word,)))
    return tuple(c for c in CATEGORIES if c in chosen)


class Geometry:
    def __init__(self, path):
        g = json.loads(path.read_text(encoding='utf-8-sig'))
        self.settings = g
        if (g['n_rsec_z'] != 1 or g['n_cry_layers'] != 1 or g['number_of_energy_bins'] != 1
                or any(v for k,v in g.items() if k.endswith('_gap') or 'virtual_crystal' in k)):
            raise ValueError('Expected one rsector ring, DOI layer and energy bin, with no gaps or virtual crystals')
        self.xy = tuple(g[f'n_{k}_xy'] for k in LEVELS)
        self.z = tuple(g[f'n_{k}_z'] for k in LEVELS)
        self.modules = g['n_rsec_xy']
        if any(type(i) is not int or i<1 for i in (self.modules,*self.xy,*self.z)):
            raise ValueError('Invalid hierarchy sizes')
        self.elements = int(np.prod(self.xy)*np.prod(self.z))
        self.detectors = self.modules*self.elements
        if g['n_det'] != self.modules*np.prod(self.xy) or g['n_rings'] != np.prod(self.z):
            raise ValueError('Hierarchy and scanner dimensions differ')
        self.energy = np.asarray([g['energy_LLD'],g['energy_ULD']],dtype=np.float32)
        self.tof = ((np.arange(g['number_of_TOF_bins']+1)-g['number_of_TOF_bins']/2)
                    *g['TOF_bin_width_mm']).astype(np.float32)
        self.block_ms = g['LM_time_block_duration_ms']
        if (type(self.block_ms) is not int or self.block_ms<1 or len(self.tof)<2
                or not np.isfinite(self.tof).all() or np.any(np.diff(self.tof)<=0)
                or not np.isfinite(self.energy).all() or np.any(np.diff(self.energy)<=0)):
            raise ValueError('Invalid time, TOF or energy binning')

    def header(self, delayeds):
        g = self.settings
        dx, dy, dz = (float(g[f'detector_{k}_dim']) for k in 'xyz')
        doi = float(g.get('average_DOI',7.))
        radius = float(g['radius'])+doi
        if not np.isfinite([dx,dy,dz,radius,doi]).all() or min(dx,dy,dz,radius)<=0 or not 0<=doi<=dx:
            raise ValueError('Invalid crystal dimensions, radius or average DOI')
        box = construct_cuboid_box((dx,dy,dz),(dx/2-doi,0.,0.),1)
        transforms=[]
        nx=int(np.prod(self.xy))
        nz=int(np.prod(self.z))
        axes=[range(n) for z,x in zip(self.z,self.xy) for n in (z,x)]
        for mz,mx,sz,sx,cz,cx in product(*axes):
            x=(mx*self.xy[1]+sx)*self.xy[2]+cx
            z=(mz*self.z[1]+sz)*self.z[2]+cz
            matrix=np.eye(4,dtype=np.float32)[:3]
            matrix[1,3]=(x-(nx-1)/2)*dy
            matrix[2,3]=(z-(nz-1)/2)*dz
            transforms.append(p.RigidTransformation(matrix=matrix))
        copies=[]
        for angle in np.arange(self.modules)*2*np.pi/self.modules:
            c,s=np.cos(angle),np.sin(angle)
            matrix=np.array([[c,-s,0,radius*c],[s,c,0,radius*s],[0,0,1,0]],dtype=np.float32)
            copies.append(p.RigidTransformation(matrix=matrix))
        elements=p.ReplicatedBoxSolidVolume(object=box,transforms=transforms)
        module=p.ReplicatedDetectorModule(object=p.DetectorModule(detecting_elements=elements),transforms=copies)
        return p.Header(scanner=p.ScannerInformation(
            model_name=g['model_name'],scanner_geometry=p.ScannerGeometry(replicated_modules=[module]),
            bulk_materials=[p.BulkMaterial(id=1,name='GATE scintillator')],
            tof_bin_edges=[[p.BinEdges(edges=self.tof)]],tof_resolution=[[g['TOF_resolution_mm']]],
            event_energy_bin_edges=[p.BinEdges(edges=self.energy)],energy_resolution_at_511=[g['EnergyResolutionAt511']],
            prompt_event_policy=p.CoincidencePolicy.OTHER,
            delayed_event_policy=p.CoincidencePolicy.OTHER if delayeds else p.CoincidencePolicy.NONE))

    def detector_id(self, rows, side):
        fields=[rows[key+side] for key in ('rsectorID','moduleID','submoduleID','crystalID','layerID')]
        sizes=[self.modules,*(x*z for x,z in zip(self.xy,self.z)),1]
        if np.any(rows['gantryID'+side]!=0) or any(np.any((v<0)|(v>=n)) for v,n in zip(fields,sizes)):
            raise ValueError('ROOT hierarchy index outside scanner JSON')
        ids=np.zeros(len(fields[0]),dtype=np.int64)
        for field,size in zip(fields,sizes):
            ids=ids*size+field
        return ids

    def castor_to_petsird(self):
        """Map CASToR ring-major IDs to the same GATE copy order used by events."""
        ring,position=np.divmod(np.arange(self.detectors),self.settings['n_det'])
        sector,position=np.divmod(position,int(np.prod(self.xy)))
        ids=sector
        for ix,iz in zip(range(len(self.xy)),range(len(self.z))):
            x=(position//int(np.prod(self.xy[ix+1:])))%self.xy[ix]
            z=(ring//int(np.prod(self.z[iz+1:])))%self.z[iz]
            ids=ids*(self.xy[ix]*self.z[iz])+z*self.xy[ix]+x
        return ids


def normalization(header, geometry, path):
    fields={}
    for line in path.with_suffix('.Cdh').read_text(encoding='utf-8-sig').splitlines():
        key,sep,value=line.partition(':')
        if sep:
            fields[key.strip().lower()]=value.strip()
    count=int(fields['number of events'])
    if (fields['data filename']!=path.name or fields['data mode'].lower()!='normalization'
            or fields['normalization correction flag']!='1' or fields['tof information flag']!='0'
            or path.stat().st_size!=count*NCF.itemsize):
        raise ValueError('Expected a complete CASToR normalization Cdf/Cdh (float32 NCF, uint32 ID pair)')
    remap=geometry.castor_to_petsird()
    n=geometry.elements
    shape=(geometry.modules-1,n,n)
    sums=np.zeros(shape,dtype=np.float64)
    hits=np.zeros(shape,dtype=np.uint32)
    rows=np.memmap(path,dtype=NCF,mode='r',shape=(count,))
    print(f'Pooling {count:,} CASToR NCF records into rotational SGIDs',flush=True)
    for start in range(0,count,1_000_000):
        batch=rows[start:start+1_000_000]
        a,b=batch['det1'],batch['det2']
        values=batch['factor']
        if np.any(a>=len(remap)) or np.any(b>=len(remap)) or not np.isfinite(values).all() or np.any(values<0):
            raise ValueError(f'Invalid CASToR normalization record at or after {start}')
        a,b=remap[a],remap[b]
        high,low=np.maximum(a,b),np.minimum(a,b)
        delta=high//n-low//n
        keep=delta>0
        keys=(delta[keep]-1,high[keep]%n,low[keep]%n)
        efficiency=np.divide(1.,values,out=np.zeros(len(batch),dtype=np.float64),where=values>0,dtype=np.float64)
        np.add.at(sums,keys,efficiency[keep])
        np.add.at(hits,keys,1)
        if start and start%100_000_000==0:
            print(f'  {start:,} normalization records read',flush=True)
    matrices=[]
    for sgid,(total,frequency) in enumerate(zip(sums,hits)):
        values=np.divide(total,frequency,out=np.zeros_like(total),where=frequency>0).astype(np.float32)
        if not np.isfinite(values).all():
            raise ValueError('Pooled PETSIRD efficiency exceeds float32 range')
        matrices.append(p.ModulePairEfficiencies(sgid=sgid,values=values))
    header.scanner.detection_efficiencies=p.DetectionEfficiencies(
        calibration_factor=float(fields.get('calibration factor',1.)),
        method_description='Inverse CASToR NCFs averaged by angular sector difference and GATE local detector IDs',
        module_pair_sgidlut=[[[[i-j-1 for j in range(i+1)] for i in range(geometry.modules)]]],
        module_pair_efficiencies_vectors=[[matrices]])
    return {'normalization_records':count,'sgid_matrices':len(matrices),
            'populated_sgid_cells':int(np.count_nonzero(hits))}


def event_blocks(tree, geometry, *, offset_ms, duration_ms, delayed_offset_ps=None, categories=CATEGORIES, limit=None):
    fields=[key+side for side in ('1','2') for key in
            ('rsectorID','moduleID','submoduleID','crystalID','layerID','gantryID','energy','time')]
    select=tuple(categories)!=CATEGORIES
    if select:
        fields.extend((key+side) for side in ('1','2') for key in ('eventID','comptonPhantom'))
    last=-1
    pending=[]
    for rows in tree.iterate(fields,entry_stop=limit,step_size=500_000,library='np'):
        time=rows['time1']*1000.
        if not np.isfinite(time).all() or np.any(time<0) or (duration_ms is not None and np.any(time>=duration_ms)):
            raise ValueError('ROOT event outside the selected per-file timestamp range')
        ticks=(np.floor((time+offset_ms)/geometry.block_ms)*geometry.block_ms).astype(np.int64)
        if np.any(np.diff(ticks)<0) or (len(ticks) and ticks[0]<last):
            raise ValueError('ROOT events are not sorted by millisecond block')
        a,b=(geometry.detector_id(rows,side) for side in ('1','2'))
        dt=(rows['time1']-rows['time2'])*1e12
        if delayed_offset_ps is not None:
            dt+=delayed_offset_ps
        swap=a<b
        tof=np.searchsorted(geometry.tof,np.where(swap,-dt,dt)*C_MM_PER_PS/2,side='right')-1
        in_window=np.ones(len(a),dtype=bool)
        for side in ('1','2'):
            energy=rows['energy'+side].astype(np.float64)*1000
            in_window&=(energy>=geometry.energy[0])&(energy<geometry.energy[-1])
        in_window&=(tof>=0)&(tof<len(geometry.tof)-1)&(a!=b)
        if select:
            comptons=rows['comptonPhantom1'].astype(np.int64)+rows['comptonPhantom2']
            category=np.where(rows['eventID1']!=rows['eventID2'],3,np.minimum(comptons,2))  # index in CATEGORIES
            in_window&=np.isin(category,[CATEGORIES.index(c) for c in categories])
        ticks=ticks[in_window]
        high=np.maximum(a,b)[in_window]
        low=np.minimum(a,b)[in_window]
        tof=tof[in_window]
        cuts=np.r_[0,np.flatnonzero(np.diff(ticks))+1,len(ticks)]
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            if lo==hi:
                continue
            tag=int(ticks[lo])
            if last not in (-1,tag):
                yield last,pending
                pending=[]
            last=tag
            pending.extend(p.CoincidenceEvent(detection_bins=[h,l],tof_idx=k)
                           for h,l,k in zip(high[lo:hi].tolist(),low[lo:hi].tolist(),tof[lo:hi].tolist()))
    if pending:
        yield last,pending


def verify(path, geometry, report):
    counts=np.zeros(2,dtype=np.int64)
    zero_prompts=0
    last=0
    with p.BinaryPETSIRDReader(str(path)) as reader:
        scanner=reader.read_header().scanner
        replicated=scanner.scanner_geometry.replicated_modules
        if (scanner.model_name!=report['scanner'] or len(replicated)!=1
                or len(replicated[0].transforms)!=geometry.modules
                or len(replicated[0].object.detecting_elements.transforms)!=geometry.elements):
            raise ValueError('PETSIRD round-trip geometry mismatch')
        efficiency=scanner.detection_efficiencies
        matrices=None
        if report['sgid_matrices']:
            groups=efficiency.module_pair_efficiencies_vectors[0][0]
            if len(groups)!=report['sgid_matrices'] or any(group.sgid!=i for i,group in enumerate(groups)):
                raise ValueError('PETSIRD round-trip SGID mismatch')
            matrices=np.asarray([group.values for group in groups],dtype=np.float32)
            if (matrices.shape!=(geometry.modules-1,geometry.elements,geometry.elements)
                    or not np.isfinite(matrices).all() or np.any(matrices<0)):
                raise ValueError('Invalid round-trip geometric normalization matrices')
        for block in reader.read_time_blocks():
            value=block.value
            lo,hi=value.time_interval.start,value.time_interval.stop
            if lo<last or hi<=lo:
                raise ValueError('Unordered PETSIRD time blocks')
            last=hi
            streams=(value.prompt_events[0][0],value.delayed_events[0][0] if value.delayed_events else [])
            for i,events in enumerate(streams):
                if not events:
                    continue
                a,b,k=np.asarray([(*ev.detection_bins,ev.tof_idx) for ev in events]).T
                if np.any((b<0)|(a<=b)|(a>=geometry.detectors)|(k<0)|(k>=len(geometry.tof)-1)):
                    raise ValueError('Invalid PETSIRD coincidence detector order or bin index')
                counts[i]+=len(events)
                if i==0 and matrices is not None:
                    n=geometry.elements
                    delta=a//n-b//n
                    zero_prompts+=int(np.count_nonzero((delta==0)|(matrices[delta-1,a%n,b%n]==0)))
    if counts.tolist()!=[sum(f[k] for f in report['files']) for k in ('prompts','delayeds')]:
        raise ValueError('PETSIRD round-trip event totals differ')
    return {'sdk_roundtrip':'PASS','last_block_end_ms':last,
            'prompts_with_zero_pooled_efficiency':zero_prompts if matrices is not None else None}


def convert(args):
    major,minor,patch=(int(x) for x in version('petsird').split('.')[:3])
    if (major,minor)!=(0,13) or patch<3:
        raise ValueError('Install PETSIRD >=0.13.3,<0.14 for this converter')
    geometry=Geometry(args.scanner)
    if args.time_mode=='per-file' and (args.file_duration_ms<1 or args.file_duration_ms%geometry.block_ms):
        raise ValueError('File duration must be a positive multiple of the block duration')
    print(f"Selection: {', '.join(args.selection)}")
    # print the timing parameters, so that the log of a conversion says which time windows were used
    tof_half_range=geometry.settings['number_of_TOF_bins']*geometry.settings['TOF_bin_width_mm']/2
    print('Timing:')
    print(f'  time mode: per-file (each ROOT file restarts its clock at 0 and starts {args.file_duration_ms} ms after the previous one)'
          if args.time_mode=='per-file' else '  time mode: absolute (the ROOT files share one clock)')
    print(f'  time blocks: {geometry.block_ms} ms')
    print(f"  TOF window: {geometry.settings['number_of_TOF_bins']} bins of {geometry.settings['TOF_bin_width_mm']} mm, "
          f'from {-tof_half_range:g} to {tof_half_range:g} mm ({-2*tof_half_range/C_MM_PER_PS:g} to {2*tof_half_range/C_MM_PER_PS:g} ps in t1 - t2)')
    print(f'  delayed coincidences: tree "delay", {args.delayed_offset_ps:g} ps added to t1 - t2'
          if args.delayed_offset_ps is not None else '  delayed coincidences: not stored',flush=True)
    files=[Path(f'{args.root_prefix}{i}.root') for i in range(1,args.number_of_root_files+1)]
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    output=args.petsird_file
    partial=output.with_suffix('.partial')
    report_path=output.with_suffix('.petsird.json')
    if output.exists() or partial.exists() or report_path.exists():
        raise FileExistsError(output)
    header=geometry.header(args.delayed_offset_ps is not None)
    report={'sdk':version('petsird'),'scanner':header.scanner.model_name,'detectors':geometry.detectors,
            'time_mode':args.time_mode,'delayed_offset_ps':args.delayed_offset_ps,
            'normalization_records':0,'sgid_matrices':0,'files':[]}
    if args.normalization_file:
        report.update(normalization(header,geometry,args.normalization_file))
    output.parent.mkdir(parents=True,exist_ok=True)
    end=0
    try:
        with p.BinaryPETSIRDWriter(str(partial)) as writer:
            writer.write_header(header)
            for index,path in enumerate(files):
                with uproot.open(path) as root:
                    offset=index*args.file_duration_ms if args.time_mode=='per-file' else 0
                    opts={'offset_ms':offset,'duration_ms':args.file_duration_ms if args.time_mode=='per-file' else None,
                          'limit':args.max_events_per_tree}
                    prompts=root['Coincidences']
                    delayed=root['delay'] if args.delayed_offset_ps is not None else None
                    pf=iter(event_blocks(prompts,geometry,categories=args.selection,**opts))
                    df=iter(event_blocks(delayed,geometry,delayed_offset_ps=args.delayed_offset_ps,**opts)) if delayed is not None else iter(())
                    a,b=next(pf,None),next(df,None)
                    p_count=d_count=0
                    while a is not None or b is not None:
                        tag=min(row[0] for row in (a,b) if row is not None)
                        if tag<end or tag+geometry.block_ms>np.iinfo(np.uint32).max:
                            raise ValueError('ROOT clock overlaps across files or exceeds PETSIRD time range; check --time-mode')
                        pr,dr=[],[]
                        if a is not None and a[0]==tag:
                            pr=a[1]
                            a=next(pf,None)
                        if b is not None and b[0]==tag:
                            dr=b[1]
                            b=next(df,None)
                        end=tag+geometry.block_ms
                        block=p.EventTimeBlock(time_interval=p.TimeInterval(start=tag,stop=end),
                                              prompt_events=[[pr]],delayed_events=[[dr]] if delayed is not None else [])
                        writer.write_time_blocks([p.TimeBlock.EventTimeBlock(block)])
                        p_count+=len(pr)
                        d_count+=len(dr)
                    item={'file':str(path),'prompt_entries':min(prompts.num_entries,args.max_events_per_tree or prompts.num_entries),
                          'delayed_entries':min(delayed.num_entries,args.max_events_per_tree or delayed.num_entries) if delayed is not None else 0,
                          'prompts':p_count,'delayeds':d_count}
                    report['files'].append(item)
                    print(item,flush=True)
            writer.write_time_blocks([])
        if args.verify:
            del header
            report.update(verify(partial,geometry,report))
        partial.rename(output)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'Saved {output} and {report_path}',flush=True)


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--root-prefix',required=True)
    cli.add_argument('--number-of-root-files',type=int,required=True)
    cli.add_argument('--petsird-file',required=True,type=Path)
    cli.add_argument('-s','--scanner',required=True,type=Path)
    cli.add_argument('-c','--normalization-file',type=Path)
    cli.add_argument('--time-mode',choices=('absolute','per-file'),required=True)
    cli.add_argument('--file-duration-ms',type=int,default=10_000)
    cli.add_argument('--delayed-offset-ps',type=float)
    cli.add_argument('--selection',type=selection,default='all',
                     help='comma-separated: true, single-scatter, multiple-scatter, random, scatter, all (default)')
    cli.add_argument('--max-events-per-tree',type=int)
    cli.add_argument('--verify',action='store_true')
    args=cli.parse_args()
    if args.number_of_root_files<0 or (args.max_events_per_tree is not None and args.max_events_per_tree<1):
        cli.error('File count must be nonnegative and event limit positive')
    if args.delayed_offset_ps is not None and not np.isfinite(args.delayed_offset_ps):
        cli.error('Delayed offset must be finite')
    convert(args)


if __name__=='__main__':
    main()
