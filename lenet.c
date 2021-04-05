#include "lenet.h"
#include "model.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

/*
DSPLIB_DATA(input,4)
_iq31 input[128];

DSPLIB_DATA(weight,4)
_iq31 weight[128];

*/
DSPLIB_DATA(mac_result,4)
_iq31 mac_result[1];

DSPLIB_DATA(result,4)
_q15 result[10];


unsigned short int i,j,z,o0,o1,i0,i1,l0,l1,x,w0,w1,y,k,tl1,tl2,length1,length2;
int temp1,temp2,temp3;

msp_mac_q15_params macParams;
msp_fill_q15_params fillParams;
msp_add_q15_params addParams;
msp_iq31_to_q15_params convParams;

msp_cmplx_fft_q15_params fftParams;
msp_cmplx_mpy_q15_params mpyParams;
msp_shift_q15_params shiftparams;
msp_cmplx_q15_params comparam;

#define S             128

#define GETLENGTH(array) (sizeof(array)/sizeof(*(array)))

#define GETCOUNT(array)  (sizeof(array)/sizeof(float))

#define FOREACH(i,count) for (i = 0; i < count; ++i)

#define CONVOLUTE_VALID(input,output,weight)                                            \
{                                                                                       \
    FOREACH(o0,GETLENGTH(output))                                                       \
        FOREACH(o1,GETLENGTH(*(output)))                                                \
            FOREACH(w0,GETLENGTH(weight))                                               \
                FOREACH(w1,GETLENGTH(*(weight)))                                        \
                    (output)[o0][o1] += (input)[o0 + w0][o1 + w1] * (weight)[w0][w1];   \
}



#define CONVOLUTE_VALID_ACCE(input,output,weight)                                            \
{                                                                                       \
    int w0l = GETLENGTH(weight);                         \
    int w1l = GETLENGTH(*(weight));                       \
    dma_transfer_macro(&(weight)[0][0],&we[0],w0l*w1l);\
    FOREACH(o0,GETLENGTH(output))                                                       \
        FOREACH(o1,GETLENGTH(*(output))){                                                \
            FOREACH(w0,w0l)                                               \
                memcpy(&in[w0*w0l],&(input)[o0 + w0][o1],w1l*sizeof(int16_t)); \
            MAC_(w0l,w1l);   \
            output[o0][o1] += res; \
        }\
}

#define CONVOLUTE_VALID_ACCE_old(input,output,weight)                                            \
{                                                                                       \
    int w0l = GETLENGTH(weight);                         \
    int w1l = GETLENGTH(*(weight));                       \
    memcpy(&we[0],&(weight)[0][0],w0l*w1l*sizeof(int16_t));    \
    FOREACH(o0,GETLENGTH(output))                                                       \
        FOREACH(o1,GETLENGTH(*(output))){                                                \
            FOREACH(w0,w0l)                                               \
                memcpy(&in[w0*w0l],&(input)[o0 + w0][o1],w1l*sizeof(int16_t));    \
            MAC_(w0l,w1l);   \
            memcpy(&(output)[o0][o1],&res,sizeof(float));    \
        }\
}

#define CONVOLUTE_FULL(input,output,weight)                                             \
{                                                                                       \
    FOREACH(i0,GETLENGTH(input))                                                        \
        FOREACH(i1,GETLENGTH(*(input)))                                                 \
            FOREACH(w0,GETLENGTH(weight))                                               \
                FOREACH(w1,GETLENGTH(*(weight)))                                        \
                (output)[i0 + w0][i1 + w1] += (input)[i0][i1] * (weight)[w0][w1];   \
}

#define CONVOLUTION_FORWARD(input,output,weight,bias,action)                    \
{                                                                               \
    for (x = 0; x < GETLENGTH(weight); ++x)                                 \
        for (y = 0; y < GETLENGTH(*weight); ++y)                            \
        CONVOLUTE_VALID(input[x], output[y], weight[x][y]);                 \
    FOREACH(j, GETLENGTH(output))                                               \
        FOREACH(i, GETCOUNT(output[j]))                                         \
        ((float *)output[j])[i] = action(((float *)output[j])[i] + bias[j]);  \
}


#define CONVOLUTION_FORWARD_NBIAS(input,output,weight,bias,action)                    \
{                                                                               \
    for (x = 0; x < GETLENGTH(weight); ++x)                                 \
        for (y = 0; y < GETLENGTH(*weight); ++y)                            \
            CONVOLUTE_VALID_ACCE(input[x], output[y], weight[x][y]);                 \
    add_bias(&output[0][0][0],bias,GETLENGTH(output),GETLENGTH(**output),GETLENGTH(*output)); \
}


#define CONVOLUTION_FIRST_PART(input,output,weight,bias,action,floop,sloop,finit,sinit)                    \
{ \
    for (x = finit; x <floop ; x++)                                 \
        for (y = sinit; y <sloop ; y++)                            \
            CONVOLUTE_VALID(input[x], output[y], weight[x%2][y]);                 \
}


#define SUBSAMP_MAX_FORWARD(input,output)                                                       \
{                                                                                               \
    const int len0 = GETLENGTH(*(input)) / GETLENGTH(*(output));                                \
    const int len1 = GETLENGTH(**(input)) / GETLENGTH(**(output));                              \
    FOREACH(i, GETLENGTH(output))                                                               \
    FOREACH(o0, GETLENGTH(*(output)))                                                           \
    FOREACH(o1, GETLENGTH(**(output)))                                                          \
    {                                                                                           \
        int x0 = 0, x1 = 0, ismax;                                                              \
        FOREACH(l0, len0)                                                                       \
            FOREACH(l1, len1)                                                                   \
        {                                                                                       \
            ismax = input[i][o0*len0 + l0][o1*len1 + l1] > input[i][o0*len0 + x0][o1*len1 + x1];\
            x0 += ismax * (l0 - x0);                                                            \
            x1 += ismax * (l1 - x1);                                                            \
        }                                                                                       \
        output[i][o0][o1] = input[i][o0*len0 + x0][o1*len1 + x1];                               \
    }                                                                                           \
}

#define DOT_PRODUCT_FORWARD(input,output,weight,bias,action)                \
{                                                                           \
    for (x = 0; x < GETLENGTH(weight); ++x)                             \
        for (y = 0; y < GETLENGTH(*weight); ++y)                        \
            output[x] += input[y] * weight[x][y];   \
    FOREACH(j, GETLENGTH(bias))                                             \
        ((float *)output)[j] = action(((float *)output)[j] + bias[j]);  \
}

#define dma_transfer_macro(input, output, size) \
{   \
    __data20_write_long((uintptr_t) &DMA0SA,(uintptr_t) input); \
    __data20_write_long((uintptr_t) &DMA0DA,(uintptr_t) output);\
    DMA0SZ = size;                       \
    DMA0CTL = DMADT_5 | DMASRCINCR_3 | DMADSTINCR_3; \
    DMA0CTL |= DMAEN;                      \
    msp_benchmarkStart(MSP_BENCHMARK_BASE, 16);\
    DMA0CTL |= DMAREQ;\
    cycleCount = msp_benchmarkStop(MSP_BENCHMARK_BASE);\
}

static void test()
{
    macParams.length = 32;
        convParams.length = 2;
        msp_mac_q15(&macParams, we, in, &mac_result[0]);
        msp_iq31_to_q15(&convParams,&mac_result[0],&res);
}

static void MAC_(int16_t w0l,int16_t w1l)
{
    res = 0;
    macParams.length = 32;
    convParams.length = 2;
    msp_mac_q15(&macParams, we, in, &mac_result[0]);
    msp_iq31_to_q15(&convParams,&mac_result[0],&res);
    //length1 += 1;
    //for(i=0;i<w0l*w1l;i++)
        //res += in[i]*we[i];
}


static void reluQ_ar(_q15 *x, int size)
{
    FOREACH(i, size)
    {
        if(x[i]<0)
            x[i] = 0;
    }
}

static void add_bias_1d(int16_t *output,const int16_t *bias,int len, bool relu)
{
    addParams.length = len;
    dma_transfer_macro(output,&circularMatrixColumn[0],len);
    dma_transfer_macro(bias,&inputMatrixColumn[0],len);
    msp_add_q15(&addParams, inputMatrixColumn, circularMatrixColumn, inputMatrixColumn);
    if(relu) reluQ_ar(inputMatrixColumn,len);
    dma_transfer_macro(inputMatrixColumn,output,len);
}

static void add_bias(int16_t *output,const int16_t *bias,int len, int len1, int len2)
{
    int slength, length;

    short int counter=0;
    FOREACH(j, len)
    {
        slength = length = len1*len2;
        memcpy(&we[0], bias + j,sizeof(int16_t));
        fillParams.value = we[0];
        fillParams.length = 256;
        msp_fill_q15(&fillParams,&circularMatrixColumn[0]);
        while(true)
        {
            if(slength<=0)
            {
                counter = 0;
                output = (output + len1*len2);
                break;
            }
            if(slength > 256)
                length = 256;
            else
                length = slength;


            addParams.length = length;


            dma_transfer_macro((output+counter*256),&inputMatrixColumn[0],length);

            msp_add_q15(&addParams, inputMatrixColumn, circularMatrixColumn, inputMatrixColumn);

            reluQ_ar(inputMatrixColumn,length);

            dma_transfer_macro(&inputMatrixColumn[0],(output+counter*256),length);

            counter ++;
            slength -= 256;
        }
    }

}

static void FC_convolution_last_old(int16_t(*action)(int16_t))
{
    for (x = 0; x < 10; ++x)
    {
        for (y = 0; y < 128; ++y)
        {
            f_output[x] += f_layerLast[y]*weight5_6[x][y];
        }
        f_output[x] = action(f_output[x] + n_bias5_6[x]);
    }
}

static void FC_convolution_last(int16_t(*action)(int16_t))
{
    macParams.length = 256;
    convParams.length = 2;

    dma_transfer_macro(&f_layerLast[0],&inputMatrixColumn[0], 256);
    for (x = 0; x < 10; ++x)
    {
        dma_transfer_macro(&weight5_6[x][0],&circularMatrixColumn[0], 256);

        msp_mac_q15(&macParams, circularMatrixColumn, inputMatrixColumn, &mac_result[0]);
        msp_iq31_to_q15(&convParams,&mac_result[0],&res);
        result[x] = res;
    }

    memcpy( &f_output[0], &result[0], 10* sizeof( int16_t ) );

    add_bias_1d(f_output,n_bias5_6,GETLENGTH(n_bias5_6),false);

}

int16_t relu(int16_t x)
{
    return x*(x > 0);
}


float relugrad(float y)
{
    return y > 0;
}


static void fftmultiplication()
{
    shiftparams.length = 128;
    //shiftparams.shift = -7;
    //comparam.length = 128;
    mpyParams.length = 128;

    fftParams.length = 128;
    fftParams.bitReverse = 1;
    fftParams.twiddleTable = msp_cmplx_twiddle_table_128_q15;

    int i,j=S*2;

    for(i=S;i>0;i--)
    {
        circularMatrixColumn[j-2] = circularMatrixColumn[i-1];
        circularMatrixColumn[j-1] = 0;

        inputMatrixColumn[j-2] = inputMatrixColumn[i-1];
        inputMatrixColumn[j-1] = 0;

        j = j-2;
    }

    //fillParams.value = 0;
    //fillParams.length = 128;

    //msp_fill_q15(&fillParams,&zerocolumn[0]);

    //msp_shift_q15(&shiftparams,inputMatrixColumn,inputMatrixColumn);
    //msp_shift_q15(&shiftparams,circularMatrixColumn,circularMatrixColumn);

    //inputMatrixColumn = {1684,1566,1090,1865,2317,815,3085,5777,0,428,7620,6446,0,6029,9002,3640,0,0,2145,0,0,759,3828,0,0,2207,392,0,734,2575,0,0,757,1540,718,0,1982,1453,0,0,0,0,1955,3567,0,0,4858,4388,10305,12467,13009,9104,6056,7891,8147,6021,2556,2749,2447,3122,2909,2745,2556,3032,3422,4303,4389,5571,2104,1772,1554,6611,809,249,4152,7550,472,2218,8611,7901,5815,7834,7873,4023,4646,6666,5130,2864,0,0,2631,583,40,2052,3067,1072,0,673,2969,3171,0,71,5168,1460,0,4029,5451,0,916,4848,1564,0,607,241,1058,0,0,0,933,1430,1447,3456,1073,6728,3246,4271,3040,6813,7608,7397,7079,5016,1838,3821,3381,3665,335,376,0,913,1104,1072,0,2122,0,0,1753,1939,420,103,8182,4639,0,364,4576,1727,0,2480,1451,49,3388,3909,3987,3902,0,0,1497,1295,0,53,1727,1157,448,1898,2378,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1813,2923,2296,2341,71,0,0,3927,1405,1386,3137,7562,1498,1344,7716,7614,4181,4664,8569,10472,2016,2794,5669,7479,1817,3367,2804,2455,2298,3391,2819,2410,4406,4658,4334,7107,2614,3320,3285,10932,2668,1742,8233,10411,2671,3633,11167,7175,2631,1510,922,1356,5789,6114,4189,3609,2191,1595,3163,4340,2051,4534,3244,4856};
    //circularMatrixColumn = {-481,-6677,7618,4423,10707,9610,2747,-7194,2386,-2903,-1650,767,-1088,1478,-1780,2358,-9931,-2158,1185,-8885,-9200,-1714,-5722,-6421,377,-4595,-3299,3569,2908,1917,-1249,-3697,-5849,1547,-585,-8200,3154,288,-4801,-8855,1790,-5158,5406,6972,-6079,2491,9196,3460,1656,1966,6235,2916,-5192,-6728,-3263,94,-319,6690,2028,-4509,-2314,-7968,6110,1457,979,6882,5889,3556,-40,-4233,582,5080,-3114,-1233,-5180,-7591,2839,-983,3761,5575,6549,-950,-1083,1551,3486,-1365,-5255,8771,2136,-634,-6299,-8185,-4398,-1786,-1368,-2994,7849,388,-77,7212,1118,-245,1620,6632,2098,4185,9156,-10221,-7484,-1684,-1528,-940,3237,3315,-711,6629,11699,536,-2905,-2424,-8120,-1334,-7746,-4996,-8144,-6865,4279,-435};

    //msp_cmplx_q15(&comparam,inputMatrixColumn,zerocolumn,inputMatrixColumn);
    //msp_cmplx_q15(&comparam,circularMatrixColumn,zerocolumn,circularMatrixColumn);

    msp_cmplx_fft_fixed_q15(&fftParams, circularMatrixColumn);

    msp_cmplx_fft_fixed_q15(&fftParams, &inputMatrixColumn[0]);         // Run LEA once for cycle count

    msp_cmplx_mpy_q15(&mpyParams, inputMatrixColumn, circularMatrixColumn, inputMatrixColumn);

    msp_cmplx_ifft_fixed_q15(&fftParams, inputMatrixColumn);

    // real number extraction needed
    i=0;
    for (z = 0; z < 128; ++z){
        inputMatrixColumn[z]=inputMatrixColumn[i];
        i=i+2;
    }

    shiftparams.shift = 7;
    msp_shift_q15(&shiftparams,inputMatrixColumn,inputMatrixColumn);
}

static void FC_Convolution()
{
    while(true)
    {
        if(length1>=512)
            break;
        int round = (length1%256);
        dma_transfer_macro(&f_layer4[0+(round/16)][0][0],&circularMatrixColumn[0], 128);
        dma_transfer_macro(&weight4_5[0+length1],&inputMatrixColumn[0], 128);


        fftmultiplication();

        dma_transfer_macro(&f_layerLast[0+round],&circularMatrixColumn[0], 128);
        addParams.length = 128;
        msp_add_q15(&addParams, inputMatrixColumn, circularMatrixColumn, inputMatrixColumn);
        dma_transfer_macro(&inputMatrixColumn[0],&f_layerLast[0+round], 128);


        length1 += 128;
    }

    add_bias_1d(f_layerLast,n_bias4_5,GETLENGTH(f_layerLast),true);
}


static void forward(int16_t(*action)(int16_t))
{
    CONVOLUTION_FORWARD_NBIAS(f_input,f_layer1, n_weight0_1, n_bias0_1, action);
    SUBSAMP_MAX_FORWARD(f_layer1, f_layer2);

    CONVOLUTION_FORWARD_NBIAS(f_layer2, f_layer3, n_weight2_3, n_bias2_3, action);
    SUBSAMP_MAX_FORWARD(f_layer3, f_layer4);

    FC_Convolution();

    //DOT_PRODUCT_FORWARD(f_layerLast, f_output, weight5_6, n_bias5_6, action);
    FC_convolution_last(action);

}




static uint8 get_result(uint8 count)
{
    int16_t *output = (int16_t *)f_output;
    //const int outlen = GETCOUNT(features->output);
    uint8 result = 0;
    int16_t maxvalue = *output;
    for (i = 1; i < count; ++i)
    {
        if (output[i] > maxvalue)
        {
            maxvalue = output[i];
            result = i;
        }
    }
    return result;
}

static void normalization()
{
    float mean = 0.5;
    float std = 0.5;
    float t = 0;
    FOREACH(j, 28)
        FOREACH(k, 28)
    {
        t = ((float)f_input[0][j][k]/255 - mean)/std;
        f_input[0][j][k] =  t * pow(2,15);
    }
}

uint8 Predict(image input,uint8 count)
{
    //fftmultiplication();
    normalization();
    forward(relu);
    return get_result(count);
}


