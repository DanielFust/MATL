

COMPILER="g++"

# Optimization Flags
OPT_FLAGS=()

if [[ "$COMPILER" == "g++" ]]; then
	ECHO "Using GNU C++ compiler"
    OPT_FLAGS+="-O"
elif [[ "$COMPILER" == "icpc" ]]; then
    ECHO "Using Intel oneAPI C++ compiler"
fi


# Debug Flags
# Custom Debug Flags
DEBUG_FLAGS=(
             "-DMATL_BOUNDS_CHECK" #Bounds checking for arrays
             "-DMATL_DIM_CHECK"    #Dimension compatibility checking for Matrix operations
            )

FLAGS=(${OPT_FLAGS[@]} ${DEBUG_FLAGS[@]})
ECHO ${FLAGS[@]}


#Include Path
INCLUDE="./include"

#Source File Path
SRC="./examples/RowMatrixDemo.cpp"

#Output Binary
OUTPUT="RowMatrixDemo.out"


#Compile
$COMPILER -std=c++17 ${FLAGS[@]} -o $OUTPUT $SRC -I $INCLUDE

#Run
./$OUTPUT

#g++ -o RowMatrix_debug.out RowMatrix_debug.cpp -I ../../../include