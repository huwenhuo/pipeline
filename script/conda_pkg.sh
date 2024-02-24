conda install ` conda list | awk '{print $1}' |grep -v '^_' | grep -v '#' | xargs conda install --yes `
