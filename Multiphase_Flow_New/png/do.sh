rm -f *.mp4 ;

ffmpeg -framerate 15 -start_number 0 -i p_50x100_Lax-Wendroff_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwp.mp4
ffmpeg -framerate 15 -start_number 0 -i d_50x100_Lax-Wendroff_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwd.mp4
ffmpeg -framerate 15 -start_number 0 -i ux_50x100_Lax-Wendroff_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwux.mp4
ffmpeg -framerate 15 -start_number 0 -i uxr_50x100_Lax-Wendroff_2d_split1_IC0_div5till5_%01dSteps.png -c:v libx264 -r 60 -pix_fmt yuv420p laxwuxr.mp4
