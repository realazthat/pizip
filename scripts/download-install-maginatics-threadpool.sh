


set -exv

mkdir -p libs
cd libs

mkdir -p maginatics-threadpool
cd maginatics-threadpool

rm -rf threadpool

git clone https://github.com/maginatics/threadpool.git

cd threadpool

git checkout d802492



