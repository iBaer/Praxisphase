#define NEQS_1D					3
#define NEQS_2D					5

#define COMP_U_1              	d
#define COMP_U_2                d * ux
#define COMP_U_3              	uxr
// 2D Y
#define COMP_U_4               	d * uy
#define COMP_U_5                uyr

#define COMP_F_1              	d*ux
#define COMP_F_2                d * ux * ux + d * cclm1 * uxr * uxr + p
#define COMP_F_3              	ux * uxr + ccl12h * uxr * uxr + powcref * gcdgcm1 * pow(p, gcm1dgc) - (p / done)
//2D Y
#define COMP_F_4               	d * ux * uy + d * cclm1 * uxr * uyr
#define COMP_F_5                ux * uyr + uy * uxr + ccl12 * uxr * uyr

#define COMP_G_1              	d * uy
#define COMP_G_2                d * ux * uy + d * cclm1 * uxr * uyr
#define COMP_G_3               	ux * uyr + uy * uxr + ccl12 * uxr * uyr
#define COMP_G_4              	d * uy * uy + d * cclm1 * uyr * uyr + p
#define COMP_G_5                uy * uyr + ccl12h * uyr * uyr + powcref * gcdgcm1 * pow(p, gcm1dgc) - (p / done)
