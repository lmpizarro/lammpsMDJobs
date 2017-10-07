rm movie.avi

for file in dump.????00.jpg; do
    echo "$file" "${file/00/}"  "${file::(-6)}.jpg"
    mv "$file"  "${file::(-6)}.jpg"
done
ffmpeg -start_number 1 -i dump.%4d.jpg -r 25 movie.avi
