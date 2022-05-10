mkdir flint_install

cd flint_install

wget https://www.flintlib.org/flint-2.8.5.tar.gz

tar xf flint-2.8.5.tar.gz

mv flint-2.8.5/* .

./configure --enable-cxx

make -j2

make check

sudo make install

cd ..

rm -rf flint_install
