rm -f *.mp4 ;

ffmpeg -framerate 5 -start_number 0 -i p_50x50_Lax-Friedrich_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 30 -pix_fmt yuv420p laxp.mp4
ffmpeg -framerate 5 -start_number 0 -i d_50x50_Lax-Friedrich_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 30 -pix_fmt yuv420p laxd.mp4
ffmpeg -framerate 5 -start_number 0 -i ux_50x50_Lax-Friedrich_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 30 -pix_fmt yuv420p laxux.mp4
ffmpeg -framerate 5 -start_number 0 -i uxr_50x50_Lax-Friedrich_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 30 -pix_fmt yuv420p laxuxr.mp4
