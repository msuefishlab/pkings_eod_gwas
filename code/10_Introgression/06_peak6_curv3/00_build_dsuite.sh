## 00_build_dsuite.sh
## Builds Dsuite with the project's Dockerfile (images/DsuiteDocker: clones Dsuite and runs make; image dsuite_local),
## then a derived image (Dockerfile.dsuite_arm; image dsuite_arm) recompiled with -fsigned-char: on aarch64 the
## unmodified build hangs in option parsing (char getopt loop). Records the version in dsuite_version.txt.
## The main Dsuite analysis used 0.5 r58; Dockerfile.dsuite_arm pins commit a547f99, which reports "0.5 r58" (the only
## commits after the r58 commit 8b0bb96 add an #include). Runs with Docker on a workstation, not on the cluster.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
outdir=${root}/output_data/10_Introgression/peak6_curv3/
mkdir -p ${outdir}

docker build -t dsuite_local -f ${root}/images/DsuiteDocker ${root}/images/
docker build -t dsuite_arm -f ${root}/code/10_Introgression/06_peak6_curv3/Dockerfile.dsuite_arm \
	${root}/code/10_Introgression/06_peak6_curv3/

{
	docker run --rm dsuite_arm bash -c "cd /opt/Dsuite && git log -1 --format='commit %H %ad' && Build/Dsuite --version && echo 'built with: make CXXFLAGS=\"-std=c++11 -fsigned-char\" (image dsuite_arm)'"
} < /dev/null > ${outdir}/dsuite_version.txt
cat ${outdir}/dsuite_version.txt
