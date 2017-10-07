rm movie.avi

for file in dump.????00.jpg; do
    echo "$file" "${file/00/}"  "${file::(-6)}.jpg"
    mv "$file"  "${file::(-6)}.jpg"
done
ffmpeg -start_number 1 -i dump.%4d.jpg -r 25 -vf scale=640x480  movie.avi
ffmpeg -start_number 20 -i dump.%4d.jpg -c:v libx264 -preset ultrafast -crf 0 output.mkv
