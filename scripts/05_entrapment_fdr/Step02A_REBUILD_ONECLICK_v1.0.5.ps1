$ErrorActionPreference = "Stop"

# ============================================================================
# Step02A REBUILD — ONE CLICK v1.0.5
# Folder-safe version: helper Python is resolved from $PSScriptRoot.
# ============================================================================

if (-not $env:ISOFORM_BENCHMARK_ROOT) {
    throw "Set ISOFORM_BENCHMARK_ROOT to the Benchmark_Program folder."
}
$ProjectRoot = (Resolve-Path -LiteralPath $env:ISOFORM_BENCHMARK_ROOT).Path

$CodeRoot = Join-Path $ProjectRoot "MANUSCRIPT_REVISION_20260813\MAIN_FIGURES_REBUILD_20260815_V01\Code"

# IMPORTANT: scripts may live inside a subfolder under Code
$ScriptDir = $PSScriptRoot
$LookupPy = Join-Path $ScriptDir "Step02A_build_lookup_v104.py"

$OutDir = Join-Path $ProjectRoot "ENTRAPMENT_FDR\Step02A_entrapment_db_v104"

$FastaName = "uniprotkb_proteome_UP000005640_2026_08_04.fasta"

Write-Host "============================================================================"
Write-Host "STEP02A ENTRAPMENT DATABASE REBUILD v1.0.5"
Write-Host "============================================================================"
Write-Host ("Script folder: {0}" -f $ScriptDir)

New-Item -ItemType Directory -Path $OutDir -Force | Out-Null

# --------------------------------------------------------------------------
# 1. Locate original FASTA
# --------------------------------------------------------------------------

Write-Host "[1/6] Locating original isoform-inclusive FASTA..."

$FastaCandidates = @(
    (Join-Path $ProjectRoot $FastaName),
    (Join-Path $CodeRoot $FastaName),
    (Join-Path $HOME "Downloads\Step01_isoform_resolvability_package\$FastaName"),
    (Join-Path $HOME "Downloads\$FastaName")
)

$Fasta = $null

foreach ($P in $FastaCandidates) {
    if (Test-Path $P) {
        $Fasta = (Get-Item $P).FullName
        break
    }
}

if ($null -eq $Fasta) {

    Write-Host "  FASTA not at expected locations; searching Desktop and Downloads..."

    $SearchRoots = @(
        [Environment]::GetFolderPath("Desktop"),
        (Join-Path $HOME "Downloads")
    )

    foreach ($R in $SearchRoots) {
        if (!(Test-Path $R)) { continue }

        $Hit = Get-ChildItem `
            -Path $R `
            -Recurse `
            -File `
            -Filter $FastaName `
            -ErrorAction SilentlyContinue |
            Select-Object -First 1

        if ($null -ne $Hit) {
            $Fasta = $Hit.FullName
            break
        }
    }
}

if ($null -eq $Fasta) {
    throw "Original FASTA not found: $FastaName"
}

Write-Host ("  FASTA: {0}" -f $Fasta)

# --------------------------------------------------------------------------
# 2. Locate FDRBench JAR
# --------------------------------------------------------------------------

Write-Host "[2/6] Locating official FDRBench JAR..."

$JarSearchRoots = @(
    $ScriptDir,
    $CodeRoot,
    (Join-Path $ProjectRoot "ENTRAPMENT_FDR"),
    (Join-Path $HOME "Downloads"),
    [Environment]::GetFolderPath("Desktop")
)

$Jar = $null

foreach ($R in $JarSearchRoots) {

    if (!(Test-Path $R)) { continue }

    $Hit = Get-ChildItem `
        -Path $R `
        -Recurse `
        -File `
        -Filter "*.jar" `
        -ErrorAction SilentlyContinue |
        Where-Object {
            $_.Name -match "FDRBench|fdrbench"
        } |
        Sort-Object Length -Descending |
        Select-Object -First 1

    if ($null -ne $Hit) {
        $Jar = $Hit.FullName
        break
    }
}

if ($null -eq $Jar) {
    throw @"
FDRBench JAR was not found.

Place the official FDRBench JAR anywhere under:
  $ScriptDir
or
  $HOME\Downloads
or
  $HOME\Desktop
and rerun this script.
"@
}

Write-Host ("  JAR:   {0}" -f $Jar)

# --------------------------------------------------------------------------
# 3. Check Java + FDRBench options
# --------------------------------------------------------------------------

Write-Host "[3/6] Checking Java and FDRBench CLI..."

& java -version
if ($LASTEXITCODE -ne 0) {
    throw "Java is unavailable or failed."
}

$Help = (& java -jar $Jar -h 2>&1 | Out-String)

if ($LASTEXITCODE -ne 0) {
    throw "Could not read FDRBench help."
}

# --------------------------------------------------------------------------
# 4. Generate target + shuffled entrapment FASTA
# --------------------------------------------------------------------------

Write-Host "[4/6] Generating 1:1 protein-level shuffled entrapment FASTA..."

$OutFasta = Join-Path $OutDir "Step02_target_plus_shuffled_entrapment_r1.fasta"

if (Test-Path $OutFasta) {
    Remove-Item $OutFasta -Force
}

$ArgsFDR = @(
    "-jar", $Jar,
    "-level", "protein",
    "-db", $Fasta,
    "-o", $OutFasta,
    "-uniprot",
    "-fix_nc", "c",
    "-check"
)

if ($Help -match "(?m)^\s*-fold\b") {
    $ArgsFDR += @("-fold", "1")
}

if (($Help -match "(?m)^\s*-seed\b") -and ($Help -match "(?m)^\s*-fix_seed\b")) {
    $ArgsFDR += @("-seed", "20260812", "-fix_seed")
}

if ($Help -match "(?m)^\s*-fix_protein_nc\b") {
    $ArgsFDR += @("-fix_protein_nc", "n")
}

Write-Host ""
Write-Host "  FDRBench command:"
Write-Host ("  java " + ($ArgsFDR -join " "))
Write-Host ""

$Log = (& java @ArgsFDR 2>&1 | Out-String)
$Log | Set-Content `
    -Path (Join-Path $OutDir "FDRBench_generation.log") `
    -Encoding UTF8

if ($LASTEXITCODE -ne 0) {
    Write-Host $Log
    throw "FDRBench database generation failed."
}

if (!(Test-Path $OutFasta)) {

    $Candidate = Get-ChildItem `
        -Path $OutDir `
        -File `
        -Filter "*.fasta" `
        -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1

    if ($null -eq $Candidate) {
        throw "FDRBench returned success but no generated FASTA was found."
    }

    $OutFasta = $Candidate.FullName
}

Write-Host ("  Generated FASTA: {0}" -f $OutFasta)

# --------------------------------------------------------------------------
# 5. Build SQLite peptide-space lookup and frozen-QC reproduction
# --------------------------------------------------------------------------

Write-Host "[5/6] Building tryptic peptide-space lookup and reproducing v103 QC..."

if (!(Test-Path $LookupPy)) {
    throw "Helper Python script not found beside the PowerShell script: $LookupPy"
}

python $LookupPy `
    --fasta $OutFasta `
    --outdir $OutDir

if ($LASTEXITCODE -ne 0) {
    throw "Peptide-space lookup/QC failed. Do NOT rerun the four workflows yet."
}

# --------------------------------------------------------------------------
# 6. Final check
# --------------------------------------------------------------------------

Write-Host "[6/6] Final files..."

$Required = @(
    "Step02_target_plus_shuffled_entrapment_r1.fasta",
    "Step02_entrapment_peptide_space.sqlite",
    "Step02A_manifest.json",
    "Step02A_space_summary.tsv",
    "Step02A_reproduction_QC.tsv"
)

$Rows = foreach ($Name in $Required) {

    $P = Join-Path $OutDir $Name

    [PSCustomObject]@{
        File = $Name
        Exists = Test-Path $P
        SizeMB = if (Test-Path $P) {
            [math]::Round((Get-Item $P).Length / 1MB, 2)
        } else {
            $null
        }
    }
}

$Rows | Format-Table -AutoSize

if (($Rows | Where-Object { -not $_.Exists }).Count -gt 0) {
    throw "One or more required Step02A output files are missing."
}

Write-Host ""
Write-Host "============================================================================"
Write-Host "STEP02A COMPLETE"
Write-Host ("Frozen entrapment database folder: {0}" -f $OutDir)
Write-Host ""
Write-Host "NEXT:"
Write-Host "Use Step02_target_plus_shuffled_entrapment_r1.fasta for AP/FP/MM/MQ reruns."
Write-Host "Keep original MBR-OFF search settings and native decoy/FDR procedures."
Write-Host "============================================================================"
