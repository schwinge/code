export CXXFLAGS="-fPIC"
export CFLAGS="-fPIC"
export NCURSES_HOME=$HOME/ncurses 
export PATH=$NCURSES_HOME/bin:$PATH
export LD_LIBRARY_PATH=$NCURSES_HOME/lib:$LD_LIBRARY_PATH
export CPPFLAGS="-I$NCURSES_HOME/include" LDFLAGS="-L$NCURSES_HOME/lib"

mkdir ncurses && cd ncurses  
wget http://ftp.gnu.org/pub/gnu/ncurses/ncurses-6.1.tar.gz  
tar -xzvf ncurses-6.1.tar.gz  
cd ncurses-6.1
./configure --prefix=$HOME/ncurses --with-shared --without-debug --enable-widec  
make && make install  

wget -O zsh.tar.xz https://sourceforge.net/projects/zsh/files/latest/download
mkdir zsh && unxz zsh.tar.xz && tar -xvf zsh.tar -C zsh --strip-components 1
cd zsh
./configure --prefix=$HOME/zsh 
make  && make install
export PATH=$HOME/zsh/bin:$PATH

sh -c "$(wget https://raw.github.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O -)"