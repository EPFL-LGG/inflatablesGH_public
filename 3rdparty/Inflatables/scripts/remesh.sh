if [ "$#" -ne 2 ]; then
    echo "usage: remesh.sh in.obj out.obj"
    exit
fi

DIR=$(dirname "$0")

meshlabserver -s $DIR/remesh_filter2.mlx -i $1 -o $2
