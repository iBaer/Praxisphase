rm -f *.mkv ;
#ffmpeg  -pattern_type glob  -i 'cont_ace_no_*.png' -vf  drawtext="text='absorbed radiation':fontsize=60:fontcolor=black:x=(main_w/2-text_w/2):y=1800" -c:v libx264 -preset ultrafast  -r 60  -pix_fmt yuv420p no.mkv &&\
ffmpeg  -pattern_type glob  -i  'cont_ace_smal*.png' -vf  drawtext="text='Grenzfläche':fontsize=60:fontcolor=black:x=(main_w/2-text_w/2):y=1800" -c:v libx264 -b:v 10M -profile:v high -preset llhq -2pass 1 -pix_fmt yuv420p  -r 60 low.mkv &
ffmpeg  -pattern_type glob  -i 'cont_ace_scale*.png' -vf  drawtext="text='Massenanteil':fontsize=60:fontcolor=black:x=(main_w/2-text_w/2):y=1800" -c:v nvenc -b:v 10M -profile:v high -preset llhq -2pass 1 -pix_fmt yuv420p  -r 60 yes.mkv 
#mplayer  -fs -loop 0  *.mkv

