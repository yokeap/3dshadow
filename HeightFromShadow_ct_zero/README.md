# Main script for 3D reconstruction from shadow on object

#### main script is divided in to two method: 
* cross inverse tracing between lower edge of shadow casted and upper edge of shadow on object.
* cross inverse tracing between lower edge of shadow casted and pseudo skeleton by means of object thinning.

```bash
├── HeightFromShadow_ct_zero
    |-- main.m             // main script based on pseudo skeleton
    |-- mainFish.m         // main script based on shadow on object
    |-- ..
    |-- .. 
    |-- (tooling functions)
```