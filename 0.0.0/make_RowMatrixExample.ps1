


#Specify C++ Compiler
$CPP_COMPILER='g++ -std=c++11'

#Path(s) to source files
$SRC=@('./examples/RowMatrixDemo.cpp')

#Path(s) to include directories
$INCLUDE=@('./include')

#Destination Directory and name
$DEST=@('./examples/RowMatrixDemo.exe')

#Debug Flags
$DEBUG_FLAGS=@(
               '-DBOUNDS_CHECK', #Matrix Bounds checking
               '-DDIM_CHECK'     #Dimension compatibility checking
              )

#Optimization Flags
$OPT_FLAGS=@()

#Concatenate Flags
$FLAGS=@($OPT_FLAGS)+($DEBUG_FLAGS)


#Build and run compile command
$CMD = $CPP_COMPILER + ' -fdump-tree-original ' + ' ' + @($FLAGS) + ' -o ' + @($DEST) + ' ' + @($SRC) + ' -I ' + @($INCLUDE) 
Invoke-Expression $CMD
 





