# Simulation

1. Download and install [IBAMR](https://ibamr.github.io) on your local computer. Alternatively, [docker containers](https://hub.docker.com/repository/docker/d0ckaaa/ibamr) have also been made available by Robert Wong.
2. Navigate to `simulation` folder.
3. Similar to the [official guide](https://ibamr.github.io/linking), run `cmake -DIBAMR_ROOT=/opt/ibamr/tmp/build/v0.14.0 ./` if you used the prebuilt Docker Image.
4. `make` and run the executable using `./main3d input3d`. Run-time of this example on an Intel Cascade Lake core is approximately 1-2 hours.

# Files

- `zebrafish_2D_coords.txt` &mdash; tracked fish coordinates
- `input3d` &mdash; IBAMR parameters
- `heading_angle.txt` &mdash; fish heading angle
- `eel3d.vertex` &mdash; point cloud of an example fish

# Contact

- Thomas Darveniza darveniza.t@wustl.edu
- Robert Wong robert.wong@wustl.edu
