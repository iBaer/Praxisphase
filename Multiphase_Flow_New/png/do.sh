rm -f *.mp4 ;

ffmpeg -framerate 15 -start_number 0 -i p_100x200_Lax-Friedrich_2d_split4_IC3_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwp.mp4
ffmpeg -framerate 20 -start_number 0 -i d_100x200_Lax-Friedrich_2d_split4_IC3_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxsplit_d.mp4
ffmpeg -framerate 15 -start_number 0 -i ux_100x200_Lax-Friedrich_2d_split4_IC3_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwux.mp4
ffmpeg -framerate 15 -start_number 0 -i uxr_100x200_Lax-Friedrich_2d_split4_IC3_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwuxr.mp4
