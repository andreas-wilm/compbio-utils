#!/bin/bash

# list directories in which given find expression is not found

echo_debug() {
    echo "DEBUG $@" 1>&2
}
echo_error() {
    echo "ERROR: $@" 1>&2
}

debug=0
sym_link_opts=""
path=""
opts=""
while [ "$1" != "" ]; do
    case $1 in
    	--debug )
    	    debug=1;
    	    ;;
    	-H | -L | -P )
    	    sym_link_opts="$sym_link_opts $1"
    	    ;;
        -* )
            # First time we encounter a dash other than -HLP means the
            # rest are all expressions
            expressions="$expressions $@"
            break
            ;;
        * )
            if [ -d "$1" ]; then
                path="$path $1"
            else
                echo_error "don't know how to handle $1"
            fi
            ;;
    esac
    shift
done

depth_opts=""
if echo "$expressions" | grep -q -- ' -m[ai][xn]depth'; then
    depth_opts=$(echo "$expressions" | sed -e 's,.* \(-m[ai][xn]depth [0-9]*\).*,\1,')
    expressions=$(echo "$expressions" | sed -e 's, -m[ai][xn]depth [0-9]*,,')
fi

test -n "$path" || exit 1
test -n "$expressions" || exit 1

if [ $debug -eq 1 ]; then
    #echo "DEBUG: remaing args after parsing: \$@=$@" 
    echo_debug "path=$path" 
    echo_debug "depth_opts=$depth_opts" 
    echo_debug "expressions=$expressions" 
    echo_debug "sym_link_opts=$sym_link_opts" 
fi

# min/maxdepth only affects dir
for d in $(find $sym_link_opts $path $depth_opts -type d); do
    echo_debug "checking d=$d" 
    if [ $(find $sym_link_opts $d $expressions | wc -l) -eq 0 ]; then
        echo $d;
    fi
done
