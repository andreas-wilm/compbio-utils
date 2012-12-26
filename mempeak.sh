#!/usr/bin/env bash

# Taken from http://stackoverflow.com/questions/774556/peak-memory-usage-of-a-linux-unix-process
# Alternative is valgrind --tool=massif

"$@" & # Run the given command line in the background.
pid=$! peak=0
while true; do
	sleep 1
	sample="$(ps -o rss= $pid 2> /dev/null)" || break
	let peak='sample > peak ? sample : peak'
done
echo "Memory Peak: $peak" 1>&2
