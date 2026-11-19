param(
    [string]$Root = 'F:\DATA\WCMC\ETSIhackathon',
    [Parameter(Mandatory=$true)][string]$Destination,
    [ValidateSet('auto','gpu','cpu','none')][string]$Execution = 'auto',
    [int]$Workers = 4
)
$ErrorActionPreference = 'Stop'
$Repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Path $Destination -Force | Out-Null
$Cases = @(
    @{ Name='IEC_6min'; Prefix='ROOT_files\IEC\ETSIPETscanner_mIEC_';
       Geometry='ETSIPETscanner_geometry.json'; Norm='Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf'; End=360 },
    @{ Name='Hoffman_3min'; Prefix='ROOT_files\voxBrain\ETSIPETscanner256mmAFOV_positronSource\ETSIPETscanner256mmAFOV_voxBrain_';
       Geometry='ETSIPETscanner256mmAFOV_geometry.json'; Norm='Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf'; End=180 },
    @{ Name='Hoffman_6min'; Prefix='ROOT_files\voxBrain\ETSIPETscanner256mmAFOV_positronSource\ETSIPETscanner256mmAFOV_voxBrain_';
       Geometry='ETSIPETscanner256mmAFOV_geometry.json'; Norm='Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf'; End=360 }
)
foreach ($Case in $Cases) {
    $Arguments = @(
        (Join-Path $Repo 'python\root_to_petsird.py'),
        '--root-prefix', (Join-Path $Root $Case.Prefix), '--number-of-root-files', '36',
        '-s', (Join-Path (Join-Path $Root 'GEOM_DATA') $Case.Geometry),
        '-c', (Join-Path (Join-Path $Root 'NORM_files') $Case.Norm),
        '--petsird-file', (Join-Path $Destination ($Case.Name+'_prompts_delayeds_norm.petsird')),
        '--execution', $Execution, '--workers', "$Workers", '--time-mode', 'auto',
        '--file-duration-ms', '10000', '--selection', 'all', '--delayeds', 'required',
        '--delayed-offset-ps', '500000', '--start-time-s', '0', '--end-time-s', "$($Case.End)"
    )
    python @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Conversion failed for $($Case.Name)" }
}
Write-Output "Completed three Python conversions in $Destination"
