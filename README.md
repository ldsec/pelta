# LattiCom

### Setting ring parameters
The ring parameters for the experiment can be adjusted in `relations/config.go`.
`mainRing` denotes the FHE ring parameters, and `rebaseRing` denotes the commitment ring parameters. 
For `rebaseRing`, we usually use the same configuration as `mainRing` -- with only `logD` parameter changed.

Note that the files need to be regenerated for every ring. Please remove the cached files in the root folder 
by invoking `rm *.test` after switching the parameters.

We give below some example configurations.

#### Zero-level PN11 with Commt. Ring Size 128
```go
var mainRing fastmath.RingParams = fastmath.BFVZeroLevelRingCustom(bfv.PN11QP54, 11)
var rebaseRing fastmath.RingParams = fastmath.BFVZeroLevelRingCustom(bfv.PN11QP54, 7)
```

#### Zero-level PN13 with Commt. Ring Size 128
```go
var mainRing fastmath.RingParams = fastmath.BFVZeroLevelRingCustom(bfv.PN13QP218, 13)
var rebaseRing fastmath.RingParams = fastmath.BFVZeroLevelRingCustom(bfv.PN13QP218, 7)
```

#### Full PN13 with Commt. Ring Size 128
```go
var mainRing fastmath.RingParams = fastmath.BFVFullRingCustom(bfv.PN13QP218, 13)
var rebaseRing fastmath.RingParams = fastmath.BFVFullRingCustom(bfv.PN13QP218, 7)
```

#### Zero-level PN15 with Commt. Ring Size 128
```go
var mainRing fastmath.RingParams = fastmath.BFVZeroLevelRingCustom(bfv.PN15QP880, 15)
var rebaseRing fastmath.RingParams = fastmath.BFVZeroLevelRingCustom(bfv.PN15QP880, 7)
```
