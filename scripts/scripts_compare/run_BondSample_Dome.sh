#!/bin/bash
#module load python3/recommended
#module load python3/recommended
get_latest_restart() {
  # $1 = restart directory
  local d="$1"

  # Prefer largest "time" in file name run_T<time>.restart
  local by_name=""
  by_name=$(
    shopt -s nullglob
    for f in "$d"/run_T*.restart; do
      # strip directory
      local base=${f##*/}
      # Match run_T<digits>[.<digits>].restart
      if [[ $base =~ ^run_T([0-9]+)(\.([0-9]+))?\.restart$ ]]; then
        local int=${BASH_REMATCH[1]}
        local frac=${BASH_REMATCH[3]}   # may be empty
        # Build key by concatenating int+frac (dot removed)
        local key="${int}${frac}"
        printf '%s %s\n' "$key" "$f"
      fi
    done | LC_ALL=C sort -k1,1n | awk '{print $2}' | tail -n 1
  )

  if [[ -n "${by_name}" ]]; then
    printf '%s\n' "${by_name}"
    return 0
  fi

  # Fallback: newest by modification time
  local by_mtime=""
  if command -v find >/dev/null 2>&1; then
    by_mtime=$(find "$d" -maxdepth 1 -type f -name 'run_T*.restart' -printf '%T@ %p\n' 2>/dev/null \
      | LC_ALL=C sort -nr \
      | awk '{ $1=""; sub(/^ /,""); print }' \
      | head -n 1 || true)
  fi

  if [[ -n "${by_mtime}" ]]; then
    printf '%s\n' "${by_mtime}"
    return 0
  fi

  # Nothing found
  return 1
}
Time_Eq=1e6
Time_Str=0
Time_Rel=0 #3e6
NumFrames=10
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=5000
tstep=0.001
Xstretch=51
MD7s=0.65
KangNC1=4.0
kappa_den=0.01
InitialVol=59167
MPs=(0.09) # 0.01 0.03 0.06 0.09 0.18 0.27 0.6  0.81
NumMons=(10 10 10 10 10 10 10 10) # 6 9 10 15 23 30 41 52
input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/RiccyProject/Revision2/Assemble_Dome/input/data_dome'
#input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/RiccyProject/Revision2/Assemble_Dome/input/data_dome_nomols'
for seed in 2 #2 3 #2 3
do
for BD in 1.35
do
for BD7s in 0.95
do
for FactorMult in 0.66
do
for kappa_den in 0.01 #5.0
do
for OvaR in 3.5 
do
for i in "${!MPs[@]}"; do
MP="${MPs[$i]}"
N_preferred="${NumMons[$i]}"
MP7s=${MP}
base='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/RiccyProject/Revision2/Assemble_Dome/runs_GCE/run_kappaD0.01_teq3e6_Frames100_tstep0.001_Nev500_KangNC14.0_MonsPref'${N_preferred}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}
restart_dir="${base}/restart"
# Pick the restart file automatically
if ! restart_path="$(get_latest_restart "${restart_dir}")"; then
echo "ERROR: No restart files found in ${restart_dir}" >&2
continue
fi
input="${restart_path}"   # full path to the chosen restart
echo "${restart_path}"
seed=1
foldernameadd='runs_BondSample/run_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MonsPref'${N_preferred}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}
#foldernameadd='runs_test/runT2_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MonsPref'${N_preferred}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
python3 build_BondSample.py ${foldernameadd} ${Time_Eq} ${Time_Str} ${Time_Rel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MD7s} ${BD7s} ${MP} ${MP7s} ${FactorMult} ${KangNC1} ${N_preferred} ${Xstretch} ${OvaR} ${kappa_den} ${InitialVol} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done
done
done
done
done