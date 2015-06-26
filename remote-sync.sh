#!/bin/sh

##USER=$1
##HOST=$2
USER="pazouki"
HOST="euler"
if [[ $USER == "" || $HOST == "" ]] ; then
	echo "usage: `basename $0` user host"
	exit 1
fi
shift 2

LOCAL_DIR=~/Repos/Git/smarticles
REMOTE_DIR=/home/pazouki/chronoStuff/smarticles

# Generic way
# REMOTE_DIR=`echo $LOCAL_DIR | sed -e "s/Users/home/"`

SAFE_DIR=/home/arman/Repos/Git/smarticles
if [[ $LOCAL_DIR/* != $SAFE_DIR/* ]] ; then
	echo "`basename $0`: not under $SAFE_DIR, won't sync"
	exit 1
fi

rsync --verbose --recursive --times \
	--exclude ".*" --exclude "cmake/" --exclude "CMakeLists.txt" \
		$LOCAL_DIR/ $USER@$HOST:$REMOTE_DIR
