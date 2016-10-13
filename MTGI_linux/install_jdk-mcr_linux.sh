mkdir $HOME
cp  jdk-8u102-linux-x64.tar.gz  $HOME
cd $HOME 
tar -zxvf jdk-8u102-linux-x64.tar.gz
rm -rf jdk-8u102-linux-x64.tar.gz
cd -
unzip MCRInstaller.zip -d MCRInstaller
cd MCRInstaller
./install -mode silent -agreeToLicense yes -destinationFolder $HOME/matlab
