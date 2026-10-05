#!/bin/sh
# Regression test driver for atompaw (used by ctest, see tests/CMakeLists.txt).
#
# Example mode: run an example input and compare the energies in the summary
# file it writes against the summary file shipped with the example.
#   run_example.sh EXE WORKDIR INPUT SUMMARY REFERENCE
#
# Value mode: run an input and compare the all-electron total energy (Ry)
# with a reference value, e.g. from the NIST atomic reference data.
#   run_example.sh EXE WORKDIR INPUT SUMMARY =VALUE TOLERANCE
#
# The all-electron total energy is the first "Total ... :" line of the summary
# file; the PAW energies are the "Summary of PAW energies" lines. Both the
# current and the older summary layouts are understood.

set -u
exe=$1; workdir=$2; input=$3; summary=$4; ref=$5

# Tolerances in Ry. The shipped example outputs were generated with v4.0.0.0,
# against which the all-electron energies agree to ~1e-8 Ry and the PAW
# energies to ~2e-6 Ry.
ae_tol=1e-7
paw_tol=1e-5

rm -rf "$workdir"
mkdir -p "$workdir" || exit 1
cp "$input" "$workdir/" || exit 1
cd "$workdir" || exit 1

"$exe" < "$(basename "$input")" > out 2> err
status=$?
if [ $status -ne 0 ]; then
  echo "FAIL: $exe exited with status $status; see $workdir/out and $workdir/err"
  tail -5 err
  exit 1
fi
if [ ! -f "$summary" ]; then
  echo "FAIL: no summary file $summary was written"
  exit 1
fi

# First all-electron total energy in a summary file
ae_total() {
  awk '$1 == "Total" && $2 == ":" { print $3; exit }' "$1"
}

# A PAW energy term: the last field of the line naming it
paw_term() {
  awk -v key="$2" 'index($0, key) { print $NF; exit }' "$1"
}

# compare NAME VALUE REFERENCE TOLERANCE -> prints a line, returns 1 on failure
compare() {
  awk -v name="$1" -v x="$2" -v r="$3" -v tol="$4" 'BEGIN {
    if (x == "" || r == "") { printf "  %-22s missing (got \"%s\", reference \"%s\")\n", name, x, r; exit 1 }
    d = x - r; if (d < 0) d = -d
    printf "  %-22s %22.12f %22.12f %10.2e %s\n", name, x, r, d, (d <= tol ? "ok" : "FAIL")
    exit (d <= tol ? 0 : 1)
  }'
}

printf "  %-22s %22s %22s %10s\n" "quantity" "this build" "reference" "|diff|"
fail=0
case "$ref" in
  =*)
    compare "AE total energy" "$(ae_total "$summary")" "${ref#=}" "$6" || fail=1
    ;;
  *)
    compare "AE total energy" "$(ae_total "$summary")" "$(ae_total "$ref")" $ae_tol || fail=1
    for key in "Total valence energy" "Smooth energy" "One center" "Smooth kinetic" \
               "Vloc energy" "Smooth exch-corr" "One-center xc"; do
      compare "$key" "$(paw_term "$summary" "$key")" "$(paw_term "$ref" "$key")" $paw_tol || fail=1
    done
    ;;
esac
exit $fail
