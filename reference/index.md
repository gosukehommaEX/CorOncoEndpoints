# Package index

## Treatment groups

Define a treatment group and calibrate the generator

- [`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
  : Define a treatment group for correlated PFS, OS and response
- [`print(`*`<OncoArm>`*`)`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/print.OncoArm.md)
  : Print an OncoArm object

## Data generation

Generate patient-level data and data at analysis cutoffs

- [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md)
  : Generate correlated PFS, OS and response
- [`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
  : Data observed at an analysis cutoff
- [`EventTime()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/EventTime.md)
  : Calendar time of a given number of events

## Properties of a group

Implied correlations, survival functions and quantiles

- [`CorEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorEndpoints.md)
  : Pearson correlations among PFS, OS and response
- [`SurvEndpoint()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/SurvEndpoint.md)
  : Survival, density and hazard of PFS or OS
- [`QuantileEndpoint()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/QuantileEndpoint.md)
  : Quantiles of PFS or OS
- [`CorBoundPFSResponse()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorBoundPFSResponse.md)
  : Attainable range of the correlation between PFS and response

## Trial design

Expected events, average hazard ratio and required events

- [`ExpectedEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md)
  : Expected number of events by calendar time
- [`AverageHR()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md)
  : Average hazard ratio between two groups
- [`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)
  : Number of events required by the log-rank test
