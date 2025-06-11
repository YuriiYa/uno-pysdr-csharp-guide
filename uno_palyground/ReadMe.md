# Getting Started

Welcome to the Uno Platform!

To discover how to get started with your new app: https://aka.platform.uno/get-started

For more information on how to use the Uno.Sdk or upgrade Uno Platform packages in your solution: https://aka.platform.uno/using-uno-sdk

Plotting

Using [scottplot Quick start](https://scottplot.net/quickstart/unoplatform/)

``` bash
dotnet add package ScottPlot
dotnet add package ScottPlot.WinUI
```

[FftSharp](https://github.com/swharden/FftSharp/blob/main/README.md)

```bash
dotnet add package FftSharp 
```

[Software Defined Radio with HackRF sdr](https://greatscottgadgets.com/sdr/) 


## Hackrf Info

- In a new PowerShell window, list connected USB devices:
    `usbipd list`
```shell
    Connected:
    BUSID  VID:PID    DEVICE                                                        STATE
    1-6    1d50:6089  HackRF One                                                    Not shared
    2-3    30c9:000e  HP Wide Vision HD Camera                                      Not shared
    2-10   8087:0026  Intel(R) Wireless Bluetooth(R)  
```

Show serial number:
`hackrf_info` 
