# xsection.dat
# 
# Tables of relative geometric properties for rounded cross-sections from 
# EPA SWMM.
#
#   Tables are named using the following convention:
#     A_xxx   = area / full area v. depth / full depth for shape xxx
#     Y_xxx   = depth / full depth v. area / full area for shape xxx
#     W_xxx   = width / max. width v. depth / full depth for shape xxx
#     R_xxx   = hyd. radius / hyd. radius at full depth v. depth /
#               full depth for shape xxx
#     S_xxx   = section factor / section factor at full depth v. area /
#               area at full depth for shape xxx
#     N_Z_xxx = number of equal intervals between 0.0 and 1.0 in the
#               table for parameter Z and shape xxx.
#

################################################################################
# CIRCULAR SHAPE
################################################################################

const    N_A_Circ   = 51
A_Circ =     # A/Afull v. Y/Yfull
[0.0,.00471,.0134,.024446,.0374,.05208,.0680,.08505,.1033,.12236,
     .1423,.16310,.1845,.20665,.2292,.25236,.2759,.29985,.3242,.34874,
     .3736,.39878,.4237,.44907,.4745,.500,.5255,.55093,.5763,.60135,
     .6264,.65126,.6758,.70015,.7241,.74764,.7708,.79335,.8154,.83690,
     .8576,.87764,.8967,.91495,.9320,.94792,.9626,.97555,.9866,.99516,1.000]

const    N_R_Circ   = 51
R_Circ =    # R/Rfull v. Y/Yfull
[.0100,.0528,.1048,.1556,.2052,.254,.3016,.3484,.3944,.4388,
 .4824,.5248,.5664,.6064,.6456,.6836,.7204,.7564,.7912,.8244,
 .8568,.888,.9176,.9464,.9736,1.00,1.024,1.048,1.070,1.0912,1.110,1.1272,
 1.144,1.1596,1.174,1.1848,1.194,1.2024,1.210,1.2148,
 1.217,1.2172,1.215,1.2104,1.203,1.192,1.178,1.1584,1.132,1.094,1.000]

const    N_Y_Circ   = 51
Y_Circ =    # Y/Yfull v. A/Afull
[0.0, 0.05236, 0.08369, 0.11025, 0.13423, 0.15643, 0.17755, 0.19772, 0.21704,
      0.23581, 0.25412, 0.27194, 0.28948, 0.30653, 0.32349, 0.34017, 0.35666,
      0.37298, 0.38915, 0.40521, 0.42117, 0.43704, 0.45284, 0.46858, 0.4843,
      0.50000, 0.51572, 0.53146, 0.54723, 0.56305, 0.57892, 0.59487, 0.61093,
      0.62710, 0.64342, 0.65991, 0.67659, 0.69350, 0.71068, 0.72816, 0.74602,
      0.76424, 0.78297, 0.80235, 0.82240, 0.84353, 0.86563, 0.88970, 0.91444,
      0.94749, 1.0]

const    N_S_Circ   = 51
S_Circ =    # S/Sfull v. A/Afull
[0.0, 0.00529, 0.01432, 0.02559, 0.03859, 0.05304, 0.06877, 0.08551, 0.10326,
      0.12195, 0.14144, 0.16162, 0.18251, 0.2041,  0.22636, 0.24918, 0.27246,
      0.29614, 0.32027, 0.34485, 0.36989, 0.39531, 0.42105, 0.44704, 0.47329,
      0.4998,  0.52658, 0.55354, 0.58064, 0.60777, 0.63499, 0.66232, 0.68995,
      0.7177,  0.74538, 0.77275, 0.79979, 0.82658, 0.8532,  0.87954, 0.90546,
      0.93095, 0.95577, 0.97976, 1.00291, 1.02443, 1.04465, 1.06135, 1.08208,
      1.07662, 1.0]

const    N_W_Circ   = 51
W_Circ =    # W/Wmax v. Y/Yfull
[0.0,   .2800,   .3919,   .4750,   .5426,   .6000,   .6499,   .6940,   .7332,
        .7684,   .8000,   .8285,   .8542,   .8773,   .8980,   .9165,   .9330, 
        .9474,   .9600,   .9708,   .9798,   .9871,   .9928,   .9968,   .9992,
        1.000,   .9992,   .9968,   .9928,   .9871,   .9798,   .9708,   .9600,
        .9474,   .9330,   .9165,   .8980,   .8773,   .8542,   .8285,   .8000,
        .7684,   .7332,   .6940,   .6499,   .6000,   .5426,   .4750,   .3919,
        .2800,   .0]

################################################################################
# EGG SHAPE
################################################################################

const    N_A_Egg   = 26
A_Egg =
[.0000,.0116,.0316,.0550,.0850,.1200,.1555,.1900,.2250,.2750,
 .3200,.3700,.4200,.4700,.5150,.5700,.6200,.6800,.7300,.7800, 
 .8350,.8850,.9250,.9550,.9800,1.000]

const    N_R_Egg   = 26
R_Egg =
[.0100,.1300,.2440,.3020,.3860,.4650,.5360,.6110,.6760,.7350, 
 .7910,.8540,.9040,.9410,1.008,1.045,1.076,1.115,1.146,1.162, 
 1.186,1.193,1.186,1.162,1.107,1.000]

const    N_Y_Egg   = 51
Y_Egg = 
[0.0, 0.04912, 0.08101, 0.11128, 0.14161, 0.16622, 0.18811, 0.21356, 0.23742,
      0.25742, 0.27742, 0.29741, 0.31742, 0.33742, 0.35747, 0.37364, 0.4,
      0.41697, 0.43372, 0.45,    0.46374, 0.47747, 0.49209, 0.50989, 0.53015,
      0.55,    0.56429, 0.57675, 0.58834, 0.6,     0.61441, 0.62967, 0.64582,
      0.66368, 0.68209, 0.7,     0.71463, 0.72807, 0.74074, 0.75296, 0.765,
      0.77784, 0.79212, 0.80945, 0.82936, 0.85,    0.86731, 0.88769, 0.914,
      0.95,    1.0]

const    N_S_Egg   = 51
S_Egg = 
[0.0, 0.00295, 0.01331, 0.02629, 0.04,    0.05657, 0.075,   0.09432, 0.11473,
      0.13657, 0.15894, 0.1803,  0.20036, 0.22,    0.23919, 0.25896, 0.28,
      0.30504, 0.33082, 0.35551, 0.37692, 0.39809, 0.42,    0.44625, 0.47321,
      0.5,     0.52255, 0.54481, 0.56785, 0.59466, 0.62485, 0.65518, 0.68181,
      0.70415, 0.72585, 0.74819, 0.77482, 0.80515, 0.83534, 0.86193, 0.88465,
      0.9069,  0.93,    0.95866, 0.98673, 1.01238, 1.03396, 1.05,    1.06517,
      1.0538,  1.0]

const    N_W_Egg   = 26
W_Egg =
[.0,    .3250,   .4270,   .5080,   .5820,   .6420,   .6960,   .7460,   .7910,
        .8360,   .8660,   .8960,   .9260,   .9560,   .9700,   .9850,   1.000,
        .9850,   .9700,   .9400,   .8960,   .8360,   .7640,   .6420,   .3100,
	.0]


################################################################################
# HORSESHOE SHAPE
################################################################################

const    N_A_Horseshoe   = 26
A_Horseshoe =
[.0000,.0181,.0508,.0908,.1326,.1757,.2201,.2655,.3118,.3587,
 .4064,.4542,.5023,.5506,.5987,.6462,.6931,.7387,.7829,.8253, 
 .8652,.9022,.9356,.9645,.9873,1.000]

const    N_R_Horseshoe   = 26
R_Horseshoe =
[.0100,.1040,.2065,.3243,.4322,.5284,.6147,.6927,.7636,.8268, 
 .8873,.9417,.9905,1.036,1.077,1.113,1.143,1.169,1.189,1.202, 
 1.208,1.206,1.195,1.170,1.126,1.000]

const    N_Y_Horseshoe   = 51
Y_Horseshoe =
[0.0, 0.04146, 0.07033, 0.09098, 0.10962, 0.12921, 0.14813, 0.16701, 0.18565,
      0.20401, 0.22211, 0.23998, 0.25769, 0.27524, 0.29265, 0.3099,  0.32704,
      0.34406, 0.36101, 0.3779,  0.39471, 0.41147, 0.42818, 0.44484, 0.46147,
      0.47807, 0.49468, 0.51134, 0.52803, 0.54474, 0.56138, 0.57804, 0.59478,
      0.61171, 0.62881, 0.64609, 0.6635,  0.68111, 0.69901, 0.71722, 0.73583,
      0.7549,  0.77447, 0.79471, 0.81564, 0.83759, 0.86067, 0.88557, 0.91159,
      0.9452,  1.0]

const    N_S_Horseshoe   = 51
S_Horseshoe =
[0.0, 0.00467, 0.01237, 0.02268, 0.03515, 0.04943, 0.06525, 0.08212, 0.10005,
      0.11891, 0.13856, 0.15896, 0.18004, 0.20172, 0.22397, 0.24677, 0.27006,
      0.2938,  0.3179,  0.34237, 0.3672,  0.39239, 0.41792, 0.44374, 0.46984,
      0.49619, 0.52276, 0.5495,  0.5764,  0.60345, 0.63065, 0.65795, 0.68531,
      0.71271, 0.74009, 0.76738, 0.79451, 0.82144, 0.84814, 0.8745,  0.90057,
      0.92652, 0.95244, 0.97724, 0.99988, 1.02048, 1.03989, 1.05698, 1.07694,
      1.07562, 1.0]

const    N_W_Horseshoe   = 26
W_Horseshoe =
[   .0,      .5878,   .8772,   .8900,   .9028,   .9156,   .9284,   .9412,
    .9540,   .9668,   .9798,   .9928,   .9992,   .9992,   .9928,   .9798,
    .9600,   .9330,   .8980,   .8542,   .8000,   .7332,   .6499,   .5426,
    .3919,   .0]

################################################################################
# GOTHIC SHAPE
################################################################################

const    N_Y_Gothic   = 51
Y_Gothic =
[0.0, 0.04522, 0.07825, 0.10646, 0.12645, 0.14645, 0.16787, 0.18641, 0.20129,
      0.22425, 0.24129, 0.25624, 0.27344, 0.29097, 0.30529, 0.32607, 0.33755,
      0.35073, 0.36447, 0.37558, 0.4,     0.4181,  0.43648, 0.45374, 0.46805,
      0.48195, 0.49626, 0.51352, 0.5319,  0.55,    0.56416, 0.57787, 0.59224,
      0.6095,  0.62941, 0.65,    0.67064, 0.69055, 0.70721, 0.72031, 0.73286,
      0.74632, 0.76432, 0.78448, 0.80421, 0.82199, 0.84363, 0.87423, 0.90617,
      0.93827, 1.0]

const    N_S_Gothic   = 51
S_Gothic =
[0.0, 0.005,   0.0174,  0.03098, 0.04272, 0.055,   0.0698,  0.0862,  0.10461,
      0.12463, 0.145,   0.16309, 0.18118, 0.2,     0.22181, 0.24487, 0.26888,
      0.2938,  0.31901, 0.34389, 0.36564, 0.38612, 0.4072,  0.43,    0.45868,
      0.48895, 0.52,    0.55032, 0.5804,  0.61,    0.63762, 0.66505, 0.6929,
      0.72342, 0.75467, 0.785,   0.81165, 0.83654, 0.86,    0.88253, 0.90414,
      0.925,   0.94486, 0.96475, 0.98567, 1.00833, 1.03,    1.0536,  1.065,
      1.055,   1.0]

const    N_W_Gothic   = 21
W_Gothic =
[0.0,   0.286,   0.643,   0.762,   0.833,   0.905,   0.952,   0.976,   0.976,
        1.0,     1.0,     0.976,   0.976,   0.952,   0.905,   0.833,   0.762,
        0.667,   0.524,   0.357,   0.0]


################################################################################
# CATENARY SHAPE
################################################################################

const    N_Y_Catenary   = 51
Y_Catenary =
[0.0, 0.02974, 0.06439, 0.08433, 0.10549, 0.12064, 0.13952, 0.1556,  0.17032,
      0.18512, 0.20057, 0.21995, 0.24011, 0.25892, 0.27595, 0.29214, 0.30802,
      0.32372, 0.33894, 0.35315, 0.36557, 0.37833, 0.3923,  0.4097,  0.42982,
      0.45,    0.46769, 0.48431, 0.5,     0.51466, 0.52886, 0.54292, 0.55729,
      0.57223, 0.5878,  0.60428, 0.62197, 0.64047, 0.6598,  0.67976, 0.7,
      0.71731, 0.73769, 0.76651, 0.8,     0.8209,  0.84311, 0.87978, 0.91576,
      0.95,    1.0]

const    N_S_Catenary   = 51
S_Catenary =
[0.0, 0.00605, 0.01455, 0.0254,  0.03863, 0.0543,  0.07127, 0.08778, 0.10372,
      0.12081, 0.14082, 0.16375, 0.18779, 0.21157, 0.23478, 0.25818, 0.28244,
      0.30741, 0.33204, 0.35505, 0.37465, 0.39404, 0.41426, 0.43804, 0.46531,
      0.49357, 0.52187, 0.54925, 0.57647, 0.60321, 0.62964, 0.65639, 0.68472,
      0.71425, 0.74303, 0.76827, 0.79168, 0.815,   0.84094, 0.86707, 0.89213,
      0.91607, 0.94,    0.96604, 0.99,    1.00714, 1.02158, 1.03814, 1.05,
      1.05,    1.0]

const    N_W_Catenary   = 21
W_Catenary =
[0.0,    0.6667,  0.8222,  0.9111,  0.9778,  1.0000,  1.0000,  0.9889,  0.9778,
         0.9556,  0.9333,  0.8889,  0.8444,  0.8000,  0.7556,  0.7000,  0.6333,
         0.5556,  0.4444,  0.3333,  0.0]


################################################################################
# SEMI-ELLIPTICAL SHAPE
################################################################################

const    N_Y_SemiEllip   = 51
Y_SemiEllip =
[0.0, 0.03075, 0.05137, 0.07032, 0.09,    0.11323, 0.13037, 0.14519, 0.15968,
      0.18459, 0.19531, 0.21354, 0.22694, 0.23947, 0.25296, 0.265,   0.27784,
      0.29212, 0.3097,  0.32982, 0.35,    0.36738, 0.3839,  0.4,     0.41667,
      0.43333, 0.45,    0.46697, 0.48372, 0.5,     0.51374, 0.52747, 0.54209,
      0.5595,  0.57941, 0.6,     0.62,    0.64,    0.66,    0.68,    0.7,
      0.71843, 0.73865, 0.76365, 0.7926,  0.82088, 0.85,    0.88341, 0.90998,
      0.93871, 1.0]

const    N_S_SemiEllip   = 51
S_SemiEllip =
[0.0, 0.00438, 0.01227, 0.02312, 0.03638, 0.05145, 0.06783, 0.085,   0.10093,
      0.11752, 0.1353,  0.15626, 0.17917, 0.20296, 0.22654, 0.24962, 0.27269,
      0.29568, 0.31848, 0.34152, 0.365,   0.38941, 0.41442, 0.44,    0.46636,
      0.49309, 0.52,    0.54628, 0.57285, 0.6,     0.62949, 0.65877, 0.68624,
      0.71017, 0.73304, 0.75578, 0.77925, 0.80368, 0.83114, 0.8595,  0.88592,
      0.90848, 0.93,    0.95292, 0.97481, 0.99374, 1.01084, 1.02858, 1.04543,
      1.05,    1.0]

const    N_W_SemiEllip   = 21
W_SemiEllip =
[0.0,    0.7000,  0.9800,  1.0000,  1.0000,  1.0000,  0.9900,  0.9800,  0.9600,
         0.9400,  0.9100,  0.8800,  0.8400,  0.8000,  0.7500,  0.7000,  0.6400,
         0.5600,  0.4600,  0.3400,  0.0]


################################################################################
# BASKETHANDLE SHAPE
################################################################################

const    N_A_Baskethandle   = 26
A_Baskethandle =
[.0000,.0173,.0457,.0828,.1271,.1765,.2270,.2775,.3280,.3780,
 .4270,.4765,.5260,.5740,.6220,.6690,.7160,.7610,.8030,.8390, 
 .8770,.9110,.9410,.9680,.9880,1.000]

const    N_R_Baskethandle   = 26
R_Baskethandle =
[.0100,.0952,.1890,.2730,.3690,.4630,.5600,.6530,.7430,.8220, 
 .8830,.9490,.9990,1.055,1.095,1.141,1.161,1.188,1.206,1.206, 
 1.206,1.205,1.196,1.168,1.127,1.000]

const    N_Y_BasketHandle   = 51
Y_BasketHandle =
[0.0, 0.04112, 0.0738,  0.1,     0.12236, 0.14141, 0.15857, 0.17462, 0.18946,
      0.20315, 0.21557, 0.22833, 0.2423,  0.25945, 0.27936, 0.3,     0.3204,
      0.34034, 0.35892, 0.37595, 0.39214, 0.40802, 0.42372, 0.43894, 0.45315,
      0.46557, 0.47833, 0.4923,  0.50945, 0.52936, 0.55,    0.57,    0.59,
      0.61023, 0.63045, 0.65,    0.66756, 0.68413, 0.7,     0.71481, 0.72984,
      0.74579, 0.76417, 0.78422, 0.80477, 0.82532, 0.85,    0.88277, 0.915,
      0.95,    1.0]

const    N_S_BasketHandle   = 51
S_BasketHandle =
[0.0, 0.00758, 0.01812, 0.03,    0.03966, 0.04957, 0.0623,  0.07849, 0.09618,
      0.11416, 0.13094, 0.14808, 0.16583, 0.18381, 0.20294, 0.225,   0.2547,
      0.28532, 0.31006, 0.32804, 0.34555, 0.36944, 0.40032, 0.43203, 0.46004,
      0.47849, 0.49591, 0.51454, 0.5381,  0.56711, 0.6,     0.64092, 0.68136,
      0.71259, 0.73438, 0.755,   0.78625, 0.8188,  0.85,    0.8679,  0.88483,
      0.90431, 0.9369,  0.97388, 1.00747, 1.033,   1.05,    1.05464, 1.06078,
      1.055,   1.0]

const    N_W_BasketHandle   = 26
W_BasketHandle =
[0.0, 0.49,    0.667,   0.82,    0.93,    1.00,    1.00,    1.00,    0.997,
      0.994,   0.988,   0.982,   0.967,   0.948,   0.928,   0.904,   0.874,
      0.842,   0.798,   0.75,    0.697,   0.637,   0.567,   0.467,   0.342,
      0.0]


################################################################################
# SEMI-CIRCULAR SHAPE
################################################################################

const    N_Y_SemiCirc   = 51
Y_SemiCirc =
[0.0, 0.04102, 0.07407, 0.1,     0.11769, 0.13037, 0.14036, 0.15,    0.16546,
      0.18213, 0.2,     0.22018, 0.2403,  0.25788, 0.27216, 0.285,   0.29704,
      0.30892, 0.32128, 0.33476, 0.35,    0.36927, 0.38963, 0.41023, 0.43045,
      0.45,    0.46769, 0.48431, 0.5,     0.51443, 0.52851, 0.54271, 0.55774,
      0.57388, 0.59101, 0.60989, 0.63005, 0.65,    0.66682, 0.68318, 0.7,
      0.71675, 0.73744, 0.76651, 0.8,     0.8209,  0.84311, 0.87978, 0.91576,
      0.95, 1.0]

const    N_S_SemiCirc   = 51
S_SemiCirc =
[0.0, 0.00757, 0.01815, 0.03,    0.0358,  0.04037, 0.04601, 0.055,   0.07475,
      0.09834, 0.125,   0.1557,  0.18588, 0.20883, 0.223,   0.23472, 0.24667,
      0.26758, 0.29346, 0.32124, 0.35,    0.3772,  0.4054,  0.43541, 0.46722,
      0.5,     0.53532, 0.56935, 0.6,     0.61544, 0.62811, 0.6417,  0.66598,
      0.7001,  0.73413, 0.76068, 0.78027, 0.8,     0.82891, 0.85964, 0.89,
      0.9127,  0.93664, 0.96677, 1.0,     1.02661, 1.04631, 1.05726, 1.06637,
      1.06,    1.0]

const    N_W_SemiCirc   = 21
W_SemiCirc =
[0.0, 0.5488,  0.8537,  1.0000,  1.0000,  0.9939,  0.9878,  0.9756,  0.9634,
      0.9451,  0.9207,  0.8902,  0.8537,  0.8171,  0.7683,  0.7073,  0.6463,
      0.5732,  0.4756,  0.3354,  0.0]

################################################################################
# SIZES FOR STANDARD ELLIPSE SHAPES
################################################################################

const NumCodesEllipse = 23

# Minor  axis (inches)
MinorAxis_Ellipse =
[ 14,19,22,24,27,29,32,34,38,43,48,53,58,63,68,72,77,82,87,92,97,106,116]

# Major axis (inches)
MajorAxis_Ellipse =
[ 23,30,34,38,42,45,49,53,60,68,76,83,91,98,106,113, 121,128,136,143,151,
  166,180]

#  Full area (sq.ft.)
Afull_Ellipse =
[ 1.80,3.30,4.10,5.10,6.30,7.40,8.80,10.20,12.90,16.60,20.50,24.80,29.50,34.60,
  40.10,46.10,52.40,59.20,66.40,74.00,82.00,99.20,118.60]

#  Full hydraulic radius (ft)
Rfull_Ellipse =
[ 0.367,0.490,0.546,0.613,0.686,0.736,0.812,0.875,0.969,1.106,1.229,1.352,
  1.475,1.598,1.721,1.845,1.967,2.091,2.215,2.340,2.461,2.707,2.968]

################################################################################
# HORIZONTAL ELLIPSE SHAPE
################################################################################

const    N_A_HorizEllipse   = 26
A_HorizEllipse =
[ .0000,.0150,.0400,.0650,.0950,.1300,.1650,.2050,.2500,.3000,
  .3550,.4150,.4800,.5200,.5850,.6450,.7000,.7500,.7950,.8350,
  .8700,.9050,.9350,.9600,.9850,1.000]

const    N_R_HorizEllipse   = 26
R_HorizEllipse =
[ .0100,.0764,.1726,.2389,.3274,.4191,.5120,.5983,.6757,.7630,
  .8326,.9114,.9702,1.030,1.091,1.146,1.185,1.225,1.257,1.274,
  1.290,1.282,1.274,1.257,1.185,1.000]

const    N_W_HorizEllipse   = 26
W_HorizEllipse =
[ .0,.3919,.5426,.6499,.7332,.8000,.8542,.8980,.9330,.9600,
  .9798,.9928,.9992,.9992,.9928,.9798,.9600,.9330,.8980,.8542,
  .8000,.7332,.6499,.5426,.3919,.0]


################################################################################
# VERTICAL ELLIPSE SHAPE
################################################################################

const    N_A_VertEllipse   = 26
A_VertEllipse =
[ .0000,.0100,.0400,.0700,.1000,.1400,.1850,.2300,.2800,.3300,
  .3800,.4300,.4800,.5200,.5700,.6200,.6700,.7200,.7700,.8150,
  .8600,.9000,.9300,.9600,.9900,1.000]

const    N_R_VertEllipse   = 26
R_VertEllipse =
[ .0100,.1250,.2436,.3536,.4474,.5484,.6366,.7155,.7768,.8396,
  .8969,.9480,.9925,1.023,1.053,1.084,1.107,1.130,1.154,1.170,
  1.177,1.177,1.170,1.162,1.122,1.000]

const    N_W_VertEllipse   = 26
W_VertEllipse =
[ .0,.3919,.5426,.6499,.7332,.8000,.8542,.8980,.9330,.9600,
  .9798,.9928,.9992,.9992,.9928,.9798,.9600,.9330,.8980,.8542,
  .8000,.7332,.6499,.5426,.3919,.0]


################################################################################
# ARCH SHAPE
################################################################################

const NumCodesArch = 102

Yfull_Arch =     # NOTE: these are in inches
[    
     # Concrete
     11,13.5,15.5,18,22.5,26.625,31.3125,36,40,45,54,62,
     72,77.5,87.125,96.875,106.5,

     # Corrugated Steel (2-2/3 x 1/2 inch Corrugation)
     13,15,18,20,24,29,33,38,43,47,52,57,

     # Corrugated Steel (3 x 1 inch Corrugation)
     31,36,41,46,51,55,59,63,67,71,75,79,83,87,91,  # 2nd value corrected

     # Structural Plate (6 x 2 inch Corrugation - Bolted Seams 
     # 19-inch Corner Radius
     55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,
     93,95,97,100,101,103,105,107,109,111,113,115,118,119,121,

     # Structural Plate (6 x 2 inch Corrugation - Bolted Seams 
     # 31-inch Corner Radius
     112,114,116,118,120,122,124,126,128,130,132,134,
     136,138,140,142,144,146,148,150,152,154,156,158]

Wmax_Arch =      # NOTE: these are in inches
[
     # Concrete
     18,22,26,28.5,36.25,43.75,51.125,58.5,65,73,88,102,
     115,122,138,154,168.75,

     # Corrugated Steel (2 2/3 x 1/2 inch corrugation)
     17,21,24,28,35,42,49,57,64,71,77,83,

     # Corrugated Steel (3 x 1 inch Corrugation)
     40,46,53,60,66,73,81,87,95,103,112,117,128,137,142,


     # Structural Plate (6 x 2 inch Corrugation - Bolted Seams
     # 19-inch Corner Radius
     73,76,81,84,87, 92,95,98,103,106,112,114,117,123,128,131,137,139,
     142,148,150,152,154,161,167,169,171,178,184,186,188,190,197,199,

     # Structural Plate (6 x 2 inch Corrugation - Bolted Seams
     # 31-inch Corner Radius
     159,162,168,170,173,179,184,187,190,195,198,204,
     206,209,215,217,223,225,231,234,236,239,245,247]

Afull_Arch =
[    1.1,1.65,2.2,2.8,4.4,6.4,8.8,11.4,14.3,17.7,25.6,34.6,
     44.5,51.7,66,81.8,99.1,1.1,1.6,2.2,2.9,4.5,6.5,8.9,11.6,14.7,18.1,
     21.9,26,7,9.4,12.3,15.6,19.3,23.2,27.4,32.1,37,42.4,48,54.2,60.5,
     67.4,74.5,22,24,26,28,31,33,35,38,40,43,46,49,52,55,58,61,64,67,
     71,74,78,81,85,89,93,97,101,105,109,113,118,122,126,131,97,102,
     105,109,114,118,123,127,132,137,142,146,151,157,161,167,172,177,
     182,188,194,200,205,211]

Rfull_Arch =
[    0.25,0.3,0.36,0.45,0.56,0.68,0.8,0.9,1.01,1.13,1.35,
     1.57,1.77,1.92,2.17,2.42,2.65,0.324,0.374,0.449,0.499,0.598,0.723,
     0.823,0.947,1.072,1.171,1.296,1.421,0.773,0.773,1.022,1.147,1.271,
     1.371,1.471,1.570,1.670,1.770,1.869,1.969,2.069,2.168,2.268,1.371,
     1.421,1.471,1.520,1.570,1.620,1.670,1.720,1.770,1.820,1.869,1.919,
     1.969,2.019,2.069,2.119,2.168,2.218,2.268,2.318,2.368,2.418,2.493,
     2.517,2.567,2.617,2.667,2.717,2.767,2.817,2.866,2.941,2.966,3.016,
     2.792,2.841,2.891,2.941,2.991,3.041,3.091,3.141,3.190,3.240,3.290,
     3.340,3.390,3.440,3.490,3.539,3.589,3.639,3.689,3.739,3.789,3.838,
     3.888,3.938]

const    N_A_Arch   = 26
A_Arch =
[ .0000,.0200,.0600,.1000,.1400,.1900,.2400,.2900,.3400,.3900,
  .4400,.4900,.5400,.5900,.6400,.6900,.7350,.7800,.8200,.8600,
  .8950,.9300,.9600,.9850,.9950,1.000]

const    N_R_Arch   = 26
R_Arch =
[ .0100,.0983,.1965,.2948,.3940,.4962,.5911,.6796,.7615,.8364,
  .9044,.9640,1.018,1.065,1.106,1.142,1.170,1.192,1.208,1.217,
  1.220,1.213,1.196,1.168,1.112,1.000]

const    N_W_Arch   = 26 
W_Arch =
[ .0,.6272,.8521,.9243,.9645,.9846,.9964,.9988,.9917,.9811,
  .9680,.9515,.9314,.9101,.8864,.8592,.8284,.7917,.7527,.7065,
  .6544,.5953,.5231,.4355,.3195,.0]