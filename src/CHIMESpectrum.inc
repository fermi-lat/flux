// File: CHIMESpectrum.inc
// blocks of numerical data included by CHIMESpectrum.cxx

namespace {
  float energies[] = {  // standard energies used by CHIME
       10.000f,     11.200f,     12.600f,     14.100f,     15.800f,  
       17.800f,     19.900f,     22.400f,     25.100f,     28.200f,  
       31.600f,     35.400f,     39.700f,     44.600f,     50.000f,  
       56.100f,     63.000f,     70.600f,     79.200f,     88.900f,  
       99.700f,    111.900f,    125.500f,    140.800f,    158.000f,  
      177.300f,    198.900f,    223.100f,    250.300f,    280.800f,  
      315.000f,    353.400f,    396.500f,    444.800f,    499.000f,  
      559.800f,    628.000f,    704.600f,    790.500f,    886.800f,  
      994.900f,   1116.100f,   1252.200f,   1404.800f,   1576.000f,  
     1768.000f,   1983.500f,   2225.300f,   2496.500f,   2800.700f,  
     3142.000f,   3525.000f,   3954.600f,   4436.600f,   4977.300f,  
     5583.900f,   6264.400f,   7027.900f,   7884.400f,   8845.300f,  
     9923.300f,  11132.700f,  12489.500f,  14011.600f,  15719.300f,  
    17635.100f,  19784.301f,  22195.500f,  24900.600f,  27935.301f,  
    31339.900f,  35159.500f,  39444.500f,  44251.801f,  49645.000f
  };

  float fluxes[] = { // differential cosmic proton flux 

1.09e-01f, 1.21e-01f, 1.36e-01f, 1.54e-01f, 1.73e-01f, 1.93e-01f, 2.16e-01f,
2.42e-01f, 2.72e-01f, 3.04e-01f, 3.40e-01f, 3.81e-01f, 4.25e-01f, 4.75e-01f,
5.30e-01f, 5.91e-01f, 6.58e-01f, 7.29e-01f, 8.07e-01f, 8.91e-01f, 9.79e-01f,
1.07e+00f, 1.17e+00f, 1.27e+00f, 1.37e+00f, 1.47e+00f, 1.56e+00f, 1.65e+00f,
1.73e+00f, 1.81e+00f, 1.86e+00f, 1.91e+00f, 1.93e+00f, 1.93e+00f, 1.91e+00f,
1.87e+00f, 1.80e+00f, 1.71e+00f, 1.60e+00f, 1.47e+00f, 1.34e+00f, 1.19e+00f,
1.05e+00f, 9.13e-01f, 7.83e-01f, 6.63e-01f, 5.55e-01f, 4.60e-01f, 3.76e-01f,
3.06e-01f, 2.45e-01f, 1.95e-01f, 1.54e-01f, 1.21e-01f, 9.39e-02f, 7.24e-02f,
5.54e-02f, 4.22e-02f, 3.19e-02f, 2.40e-02f, 1.80e-02f, 1.34e-02f, 9.91e-03f,
7.32e-03f, 5.37e-03f, 3.94e-03f, 2.88e-03f, 2.09e-03f, 1.52e-03f, 1.10e-03f,
8.00e-04f, 5.81e-04f, 4.21e-04f, 3.08e-04f, 2.24e-04f
};

//check variable to keep flux table from being modified multiple times
int fluxcheck = 0;

  // cosmic ray proton fluxes vs. latitude and longitude
  float gfluxes[] = {
    132.165f, 104.037f,  86.241f,  74.649f,  68.091f,  64.588f,  64.211f,
     66.958f,  72.777f,  83.434f, 100.197f, 127.220f, 171.465f,
    138.386f, 107.143f,  87.991f,  75.601f,  68.240f,  64.130f,  63.167f, 
     65.285f,  70.616f,  80.125f,  95.384f, 120.158f, 159.798f,
    145.207f, 110.980f,  89.893f,  76.628f,  68.471f,  63.699f,  62.074f, 
     63.680f,  68.424f,  76.749f,  90.609f, 113.109f, 150.392f,
    152.382f, 115.636f,  92.376f,  77.850f,  68.751f,  63.289f,  61.157f, 
     61.971f,  66.149f,  73.675f,  86.653f, 106.896f, 140.820f,
    160.464f, 120.481f,  95.098f,  79.158f,  69.131f,  62.995f,  60.371f,  
     60.645f,  63.946f,  70.974f,  82.630f, 101.360f, 131.610f,
    170.523f, 125.733f,  98.038f,  80.610f,  69.537f,  62.661f,  59.616f,  
     59.374f,  61.862f,  68.210f,  78.607f,  95.773f, 123.986f,
    180.763f, 131.210f, 101.200f,  82.183f,  70.099f,  62.400f,  58.872f, 
     58.080f,  60.121f,  65.427f,  74.965f,  90.647f, 116.101f,
    191.618f, 137.972f, 104.573f,  83.923f,  70.727f,  62.292f,  58.172f,  
     56.835f,  58.363f,  62.803f,  71.818f,  86.234f, 109.314f,
    205.645f, 145.471f, 108.071f,  85.779f,  71.398f,  62.276f,  57.538f, 
     55.615f,  56.679f,  60.668f,  68.670f,  81.855f, 103.462f,
    219.581f, 153.258f, 111.875f,  87.744f,  72.189f,  62.293f,  56.903f,  
     54.473f,  54.980f,  58.584f,  65.550f,  77.617f,  97.673f,
    235.025f, 161.089f, 116.709f,  89.766f,  72.996f,  62.351f,  56.324f, 
     53.286f,  53.352f,  56.511f,  62.890f,  74.303f,  92.569f,
    253.010f, 170.734f, 121.680f,  91.866f,  73.893f,  62.435f,  55.780f, 
     52.277f,  51.968f,  54.515f,  60.567f,  71.041f,  88.259f,
    270.364f, 180.745f, 126.653f,  94.065f,  74.788f,  62.520f,  55.236f, 
     51.490f,  50.760f,  52.714f,  58.295f,  67.894f,  84.118f,
    291.127f, 189.978f, 131.510f,  96.610f,  75.675f,  62.628f,  54.734f, 
     50.698f,  49.592f,  51.237f,  56.131f,  65.006f,  80.143f,
    311.591f, 200.446f, 136.174f,  99.132f,  76.548f,  62.772f,  54.237f,  
     49.870f,  48.434f,  49.774f,  53.995f,  62.748f,  77.065f,
    332.071f, 211.291f, 141.461f, 101.435f,  77.374f,  62.795f,  53.733f,  
     49.115f,  47.302f,  48.409f,  52.477f,  60.605f,  74.293f,
    353.875f, 221.539f, 146.811f, 103.608f,  78.069f,  62.798f,  53.377f,  
     48.345f,  46.225f,  47.101f,  50.986f,  58.639f,  71.726f,
    373.163f, 230.389f, 151.387f, 105.494f,  78.648f,  62.799f,  52.967f,  
     47.593f,  45.199f,  45.825f,  49.612f,  56.738f,  69.393f,
    392.107f, 239.360f, 155.262f, 106.922f,  79.087f,  62.616f,  52.511f,  
     46.828f,  44.370f,  44.748f,  48.303f,  55.086f,  67.269f,
    409.782f, 247.231f, 158.382f, 107.981f,  79.276f,  62.378f,  52.003f,  
     46.052f,  43.627f,  43.914f,  47.152f,  53.925f,  65.835f,
    423.626f, 253.312f, 160.420f, 108.597f,  79.259f,  61.989f,  51.510f, 
     45.292f,  42.899f,  43.114f,  46.070f,  52.906f,  64.660f,
    432.302f, 257.141f, 161.538f, 108.868f,  78.951f,  61.446f,  50.863f,  
     44.712f,  42.228f,  42.406f,  45.262f,  52.057f,  63.710f,
    436.541f, 258.679f, 161.524f, 108.342f,  78.405f,  60.793f,  50.184f, 
     44.143f,  41.588f,  41.776f,  44.694f,  51.395f,  63.039f,
    436.971f, 257.662f, 160.529f, 107.310f,  77.530f,  59.961f,  49.470f, 
     43.533f,  41.027f,  41.217f,  44.225f,  50.960f,  62.659f,
    432.597f, 254.165f, 158.382f, 105.664f,  76.423f,  58.992f,  48.699f,  
     42.967f,  40.494f,  40.798f,  43.903f,  50.709f,  62.587f,
    423.437f, 247.915f, 154.876f, 103.394f,  74.983f,  57.861f,  47.886f,  
     42.411f,  40.013f,  40.438f,  43.744f,  50.668f,  62.804f,
    408.802f, 239.324f, 150.217f, 100.654f,  73.265f,  56.668f,  47.010f, 
     41.857f,  39.613f,  40.234f,  43.708f,  50.866f,  63.376f,
    388.683f, 230.310f, 144.454f,  97.534f,  71.347f,  55.484f,  46.196f, 
     41.391f,  39.361f,  40.184f,  43.898f,  51.306f,  64.257f,
    369.555f, 219.723f, 138.491f,  94.557f,  69.191f,  54.395f,  45.555f, 
     40.965f,  39.178f,  40.247f,  44.175f,  52.024f,  65.491f,
    346.290f, 207.321f, 132.438f,  91.253f,  67.253f,  53.289f,  44.934f, 
     40.656f,  39.119f,  40.455f,  44.693f,  52.990f,  67.016f,
    321.949f, 194.846f, 125.590f,  87.613f,  65.360f,  52.183f,  44.357f, 
     40.373f,  39.196f,  40.853f,  45.407f,  54.186f,  68.934f,
    299.084f, 182.772f, 118.455f,  83.868f,  63.452f,  51.093f,  43.815f,  
     40.243f,  39.395f,  41.421f,  46.254f,  55.637f,  71.692f,
    274.309f, 169.934f, 112.713f,  80.415f,  61.583f,  50.058f,  43.434f, 
     40.244f,  39.807f,  42.164f,  47.544f,  57.535f,  74.866f,
    252.445f, 159.005f, 106.775f,  77.346f,  59.760f,  49.150f,  43.127f,  
     40.358f,  40.329f,  43.080f,  49.135f,  59.979f,  78.382f,
    231.039f, 147.999f, 100.979f,  74.358f,  58.076f,  48.362f,  42.943f, 
     40.627f,  41.026f,  44.190f,  50.925f,  62.717f,  82.220f,
    212.594f, 137.667f,  95.748f,  71.498f,  56.510f,  47.706f,  42.887f,  
     41.068f,  41.903f,  45.401f,  52.954f,  65.673f,  87.318f,
    194.796f, 128.842f,  91.295f,  68.836f,  55.196f,  47.256f,  43.022f,  
     41.647f,  42.890f,  47.037f,  55.141f,  69.057f,  92.867f,
    179.955f, 120.427f,  87.061f,  66.519f,  54.231f,  46.957f,  43.251f,  
     42.365f,  44.036f,  49.000f,  57.897f,  73.311f,  98.583f,
    165.635f, 113.261f,  83.113f,  64.724f,  53.459f,  46.848f,  43.637f,  
     43.218f,  45.363f,  51.113f,  61.100f,  77.782f, 106.100f,
    154.656f, 107.324f,  79.695f,  63.182f,  52.852f,  46.913f,  44.164f,  
     44.221f,  47.107f,  53.343f,  64.378f,  82.582f, 113.609f,
    143.760f, 101.975f,  77.065f,  61.910f,  52.429f,  47.159f,  44.837f, 
     45.392f,  49.024f,  55.841f,  67.937f,  88.476f, 122.646f,
    134.531f,  97.095f,  74.798f,  60.821f,  52.230f,  47.577f,  45.830f,  
     46.949f,  51.053f,  58.959f,  72.363f,  94.452f, 132.691f,
    127.305f,  93.172f,  72.878f,  60.027f,  52.147f,  48.189f,  46.967f, 
     48.619f,  53.174f,  62.066f,  76.735f, 101.157f, 143.371f,
    120.645f,  90.017f,  71.232f,  59.455f,  52.278f,  48.890f,  48.213f,  
     50.350f,  55.582f,  65.236f,  81.377f, 108.462f, 155.734f,
    114.793f,  87.291f,  69.922f,  59.113f,  52.557f,  49.704f,  49.545f, 
     52.136f,  58.317f,  68.987f,  86.859f, 115.583f, 167.324f,
    110.379f,  85.016f,  68.886f,  58.937f,  52.965f,  50.600f,  50.927f, 
     54.015f,  61.097f,  72.912f,  92.161f, 124.691f, 181.445f,
    106.875f,  83.091f,  68.132f,  59.015f,  53.584f,  51.603f,  52.359f,  
     56.349f,  63.776f,  76.681f,  97.752f, 133.365f, 194.624f,
    103.831f,  81.584f,  67.640f,  59.195f,  54.384f,  52.639f,  54.041f, 
     58.621f,  66.810f,  80.741f, 103.808f, 142.582f, 209.578f,
    101.227f,  80.381f,  67.409f,  59.525f,  55.225f,  54.002f,  55.858f, 
     60.813f,  69.998f,  85.121f, 109.469f, 152.027f, 223.747f, 
     99.110f,  79.517f,  67.308f,  59.924f,  56.202f,  55.381f,  57.616f, 
     62.934f,  73.013f,  89.185f, 115.369f, 160.747f, 237.643f,
     97.491f,  78.923f,  67.464f,  60.435f,  57.147f,  56.738f,  59.327f,  
     65.250f,  75.785f,  92.929f, 121.817f, 169.809f, 252.184f,
     96.248f,  78.592f,  67.680f,  61.016f,  58.116f,  58.096f,  60.944f,
     67.609f,  78.698f,  97.154f, 127.501f, 178.636f, 264.596f,
     95.326f,  78.490f,  68.069f,  61.632f,  59.076f,  59.331f,  62.446f, 
     69.735f,  81.677f, 100.998f, 132.503f, 186.187f, 276.359f,
     94.783f,  78.549f,  68.482f,  62.326f,  60.009f,  60.551f,  64.175f, 
     71.705f,  84.266f, 104.395f, 137.394f, 192.234f, 287.513f,
     94.544f,  78.799f,  69.002f,  63.220f,  60.918f,  61.605f,  65.739f, 
     73.406f,  86.490f, 107.250f, 142.012f, 198.897f, 296.364f,
     94.560f,  79.191f,  69.605f,  64.095f,  61.781f,  62.784f,  67.167f, 
     74.872f,  88.395f, 109.454f, 145.599f, 204.042f, 303.026f,
     94.793f,  79.641f,  70.182f,  64.885f,  62.800f,  63.978f,  68.346f, 
     76.321f,  89.868f, 111.713f, 148.341f, 207.785f, 307.346f,
     95.291f,  80.260f,  70.880f,  65.664f,  63.668f,  64.893f,  69.358f, 
     77.584f,  91.112f, 113.446f, 150.065f, 209.868f, 309.383f,
     96.024f,  80.995f,  71.523f,  66.417f,  64.488f,  65.701f,  70.111f,
     78.512f,  92.143f, 114.478f, 151.004f, 210.479f, 309.338f,
     96.906f,  81.820f,  72.196f,  67.161f,  65.182f,  66.338f,  70.675f,
     79.101f,  92.701f, 114.944f, 151.049f, 209.822f, 307.050f,  
     97.974f,  82.691f,  72.866f,  67.796f,  65.762f,  66.878f,  71.065f,
     79.420f,  92.840f, 114.603f, 150.195f, 207.733f, 302.744f,
     99.195f,  83.639f,  73.562f,  68.401f,  66.293f,  67.260f,  71.252f, 
     79.431f,  92.495f, 113.901f, 148.519f, 204.308f, 296.654f,
    100.634f,  84.684f,  74.467f,  68.995f,  66.716f,  67.467f,  71.261f,
     79.143f,  91.793f, 112.366f, 146.118f, 199.578f, 288.711f,
    102.232f,  85.770f,  75.416f,  69.540f,  67.031f,  67.580f,  71.095f,
     78.616f,  90.711f, 110.409f, 142.994f, 193.820f, 278.584f,
    104.035f,  87.015f,  76.363f,  70.043f,  67.276f,  67.535f,  70.751f, 
     77.821f,  89.301f, 107.986f, 139.124f, 187.486f, 266.864f,
    106.065f,  88.321f,  77.287f,  70.534f,  67.454f,  67.330f,  70.305f,
     76.835f,  87.874f, 105.686f, 134.438f, 181.042f, 256.438f,
    108.506f,  89.937f,  78.326f,  71.015f,  67.605f,  67.134f,  69.640f, 
     75.611f,  86.183f, 103.036f, 129.744f, 173.915f, 244.372f,
    111.692f,  91.800f,  79.398f,  71.496f,  67.699f,  66.849f,  68.952f,
     74.202f,  84.300f, 100.074f, 125.217f, 165.956f, 231.306f,
    114.999f,  93.871f,  80.585f,  72.012f,  67.755f,  66.453f,  68.158f,
     72.831f,  82.190f,  96.801f, 120.425f, 157.520f, 218.766f, 
    118.684f,  96.105f,  81.807f,  72.555f,  67.794f,  66.042f,  67.258f,
     71.512f,  79.925f,  93.341f, 115.125f, 150.184f, 206.863f,
    122.771f,  98.519f,  83.162f,  73.133f,  67.852f,  65.577f,  66.288f, 
     70.021f,  77.521f,  89.715f, 109.526f, 142.498f, 193.852f,
    127.049f, 101.158f,  84.625f,  73.824f,  67.984f,  65.083f,  65.275f, 
     68.498f,  74.972f,  86.654f, 104.877f, 134.218f, 182.489f,
    132.165f, 104.037f,  86.241f,  74.649f,  68.091f,  64.588f,  64.211f,
     66.958f,  72.777f,  83.434f, 100.197f, 127.220f, 171.465f
  };

}   // anonymous namespace 
