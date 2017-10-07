rm movie.avi

for file in dump.???00.jpg; do
    mv "$file" "${file/00/}"
done
ffmpeg -start_number 6 -i dump.%3d.jpg -r 25 movie.avi
