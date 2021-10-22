# A Low-Latency Library in FPGA Hardware for High-Frequency Trading (HFT)

## Abstract

System used nowadays have high and unpredictable latency. A solution can be found by using FPGAs and in particular a pre-built library to accelerate software developement.

## Introduction 

HFT: Rapid electronic trade in financial instruments, using financial datas coming from trading venues and rapid buy/sell operation with small margin over a large amount of options. 
HFT is more thano 70% of trades during 2010.

### Traders
Traders use UDP/IP market data stream to get price informations from trading venues. Trading decisions are taken using complex algorithms and can be quick or slow.

### Broker

Brokers take orders from clients and passes them to the exchange. Order routing is required to get the best profit over the trade and to avoid the matching fees. Faster match and larger troughput are two featureas attractive for traders that cares about latency. 

### Low Latency
HFT Traders attemp to exploit fleeting inefficiencies in the market. With an increasing number of HFT traders the ones with the best latency are the ones that can retrieve the best profit. 

Slippage events can also happen, these ones are related to orders for which the price is worse than expected. 

## Current platforms for HFT
Usually there is software implemente on commodity servers, but it is problematic because of unpredictable response times.

### Software base HFT
Software is usually designed to have minimum latency, so the main contribution is due to computer and network stack, that can be eventually bypassed using specialized low latency cards.

Measurements on Linux server estimate a half round trip o 15-20 $\mu sec$. Results of about 2.9 and 6 $\mu s$ were measured with a TCP offload engine. This estimates * does not * take into account additional latency of end-user applications. 

We can also consider Myricom's low latency NIC that allow to send orders bypassing the system and kernel.

Infiniband devices have been introducted in the market and they provides even smaller latency, but ethernet is prevalent and will be a standard for a long time.

### Custom Hardware-based HFT platforms

Industry has explored custom hardware approaches, but ASICs lack of flexibility for new protocols and GPUs are optimized for throughput and can't have low latency.

FPGAs 
