
#include <msp430.h>
#include "lenet.h"
//#include "model.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#define CS_8MHZ             0
#define CS_16MHZ            1



short int i,percent,j,predict;
uint8 left;

image *test_data={0};
uint8 *test_label={0};


void testFirstImage()
{

    left = 7;
    predict = Predict(test_data[0], 10);

    if(left == predict)
    {
        printf("Shauo! prediction milse. Ha ha!");
    }

}


void configClocks(void)
{
#if CS_8MHZ
    //Ensure that there are 0 wait states enabled
    FRCTL0 = FRCTLPW | NWAITS_0;


    // Clock System Setup
    CSCTL0_H = CSKEY_H;                     // Unlock CS registers
    CSCTL1 = DCOFSEL_0;                     // Set DCO to 1MHz

    // Set SMCLK = MCLK = DCO, ACLK = VLOCLK
    CSCTL2 = SELA__VLOCLK | SELS__DCOCLK | SELM__DCOCLK;

    // Per Errata CS12 set divider to 4 before changing frequency to
    // prevent out of spec operation from overshoot transient
    CSCTL3 = DIVA__4 | DIVS__4 | DIVM__4;   // Set all corresponding clk sources to divide by 4 for errata
    CSCTL1 = DCOFSEL_6;                     // Set DCO to 8MHz

    // Delay by ~10us to let DCO settle. 60 cycles = 20 cycles buffer + (10us / (1/4MHz))
    __delay_cycles(60);

    CSCTL3 = DIVA__1 | DIVS__1 | DIVM__1;   // Set MCLK  to 8MHz, SMCLK to 1MHz

    CSCTL4 = HFXTOFF + LFXTOFF;             // Cancel external clocks
    CSCTL5 &= ~(HFXTOFFG + LFXTOFFG);

    CSCTL0_H = 0;                           // Lock CS Registers

#endif

#if CS_16MHZ
    // Configure one FRAM waitstate as required by the device datasheet for MCLK
    // operation beyond 8MHz _before_ configuring the clock system.
    FRCTL0 = FRCTLPW | NWAITS_1;

    // Clock System Setup
    CSCTL0_H = CSKEY_H;                     // Unlock CS registers
    CSCTL1 = DCOFSEL_0;                     // Set DCO to 1MHz

    // Set SMCLK = MCLK = DCO, ACLK = VLOCLK
    CSCTL2 = SELA__VLOCLK | SELS__DCOCLK | SELM__DCOCLK;

    // Per Errata CS12 set divider to 4 before changing frequency to
    // prevent out of spec operation from overshoot transient
    CSCTL3 = DIVA__4 | DIVS__4 | DIVM__4;   // Set all corresponding clk sources to divide by 4 for errata
    CSCTL1 = DCOFSEL_4 | DCORSEL;           // Set DCO to 16MHz

    // Delay by ~10us to let DCO settle. 60 cycles = 20 cycles buffer + (10us / (1/4MHz))
    __delay_cycles(60);

    CSCTL3 = DIVA__1 | DIVS__1 | DIVM__1;   // Set MCLK to 16MHz, SMCLK to 1MHz

    CSCTL4 = HFXTOFF + LFXTOFF;             // Cancel external clocks
    CSCTL5 &= ~(HFXTOFFG + LFXTOFFG);

    CSCTL0_H = 0;                           // Lock CS Registers
#endif
}


void configGPIO(void)
{
    // Configure ports for low power
    PADIR = 0xFFFF;  PAOUT = 0;
    PBDIR = 0xFFFF;  PBOUT = 0;
    PCDIR = 0xFFFF;  PCOUT = 0;
    PDDIR = 0xFFFF;  PDOUT = 0;
    PJDIR = 0xFFFF;  PJOUT = 0;
    P1DIR = 0xFF;  P1OUT = 0;
    P2DIR = 0xFF;  P2OUT = 0;
    P3DIR = 0xFF;  P3OUT = 0;
    P4DIR = 0xFF;  P4OUT = 0;
    P5DIR = 0xFF;  P5OUT = 0;
    P6DIR = 0xFF;  P6OUT = 0;
    P7DIR = 0xFF;  P7OUT = 0;
    P8DIR = 0xFF;  P8OUT = 0;
}

int main(void)
{
    WDTCTL = WDTPW | WDTHOLD;   // stop watchdog timer

    PM5CTL0 &= ~LOCKLPM5;
    configGPIO();
    configClocks();
    testFirstImage();

    return 0;
}


