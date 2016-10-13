export JAVA_HOME=$HOME/jdk1.8.0_102
export JRE_HOME=$HOME/jdk_1.8.0_102/jre
export JAVA_BIN=$HOME/jdk_1.8.0_102/bin
export PATH=$JAVA_HOME/bin:$PATH
export CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lob/tools.jar
export JAVA_HOME JAVA_BIN PATH CLASSPATH

export MAC_HOME=$HOME/matlab
LD_LIBRARY_PATH=$HOME/matlab/v81/bin/glnxa64
XAPPLRESDIR=$HOME/matlab/v81/x11/app-defaults
LD_PATH=$HOME/matlab/v81/runtime/glnxa64
LD_LIBRARY_PATH=$LD_PATH:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
export XAPPLRESDIR

java -jar MTGI_linux.jar
