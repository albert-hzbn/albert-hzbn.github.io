class GrainGrowthSimulation {
  constructor(params) {
    // Constants
    this.Nx = params.Nx;
    this.Ny = params.Ny;
    this.ngrains = params.ngrains;
    this.dx = params.dx;
    this.dy = params.dy;
    this.dt = params.dt;
    this.kappa = params.kappa;
    this.L = params.L;
    this.tEnd = params.tEnd;
    this.outInterval = params.outInterval;

    // Initialize arrays
    this.etas = Array(this.ngrains)
      .fill()
      .map(() =>
        Array(this.Nx)
          .fill()
          .map(() => Array(this.Ny).fill(0))
      );
    this.glist = Array(this.ngrains).fill(1);
    this.lap_etas = Array(this.Nx)
      .fill()
      .map(() => Array(this.Ny).fill(0));

    this.areaFractions = [];
  }

  distance(p1, p2) {
    return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
  }

  generateRandomPoints(numPoints) {
    const points = [];
    for (let i = 0; i < numPoints; i++) {
      points.push({
        x: Math.floor(Math.random() * this.Nx),
        y: Math.floor(Math.random() * this.Ny),
      });
    }
    return points;
  }

  generatePolycrystalline(grainCenters) {
    for (let ind_x = 0; ind_x < this.Nx; ind_x++) {
      for (let ind_y = 0; ind_y < this.Ny; ind_y++) {
        let nearestGrainIndex = 0;
        let minDistance = this.distance(grainCenters[0], {
          x: ind_x,
          y: ind_y,
        });

        // Find the nearest grain for each point
        for (let i = 1; i < this.ngrains; i++) {
          const dist = this.distance(grainCenters[i], { x: ind_x, y: ind_y });
          if (dist < minDistance) {
            minDistance = dist;
            nearestGrainIndex = i;
          }
        }

        // Set etas for the nearest grain
        for (let i = 0; i < this.ngrains; i++) {
          this.etas[i][ind_x][ind_y] = i === nearestGrainIndex ? 1.0 : 0.0;
        }
      }
    }
  }

  grainGrowthStep() {
    for (let ind_grain = 0; ind_grain < this.ngrains; ind_grain++) {
      if (this.glist[ind_grain] === 1) {
        for (let ind_x = 0; ind_x < this.Nx; ind_x++) {
          for (let ind_y = 0; ind_y < this.Ny; ind_y++) {
            // Periodic boundary
            const periodicIndex = (idx, size) => ((idx % size) + size) % size;
            const iw = periodicIndex(ind_x - 1, this.Nx);
            const ie = periodicIndex(ind_x + 1, this.Nx);
            const jw = periodicIndex(ind_y - 1, this.Ny);
            const je = periodicIndex(ind_y + 1, this.Ny);

            // Calculate Laplacian
            this.lap_etas[ind_x][ind_y] =
              (this.etas[ind_grain][ie][ind_y] +
                this.etas[ind_grain][iw][ind_y] +
                this.etas[ind_grain][ind_x][je] +
                this.etas[ind_grain][ind_x][jw] -
                4.0 * this.etas[ind_grain][ind_x][ind_y]) /
              (this.dx * this.dy);

            // Free energy derivative
            const alpha = 1.0;
            const beta = 1.0;
            let sum = 0.0;

            for (let ind_grain2 = 0; ind_grain2 < this.ngrains; ind_grain2++) {
              if (ind_grain2 !== ind_grain) {
                sum += Math.pow(this.etas[ind_grain2][ind_x][ind_y], 2);
              }
            }

            const dfdetla =
              alpha *
              (2 * beta * this.etas[ind_grain][ind_x][ind_y] * sum +
                Math.pow(this.etas[ind_grain][ind_x][ind_y], 3) -
                this.etas[ind_grain][ind_x][ind_y]);

            this.etas[ind_grain][ind_x][ind_y] =
              this.etas[ind_grain][ind_x][ind_y] -
              this.dt *
                this.L *
                (dfdetla - this.kappa * this.lap_etas[ind_x][ind_y]);

            // For small deviations
            if (this.etas[ind_grain][ind_x][ind_y] > 1.0) {
              this.etas[ind_grain][ind_x][ind_y] = 1.0;
            }
            if (this.etas[ind_grain][ind_x][ind_y] < 0.0) {
              this.etas[ind_grain][ind_x][ind_y] = 0.0;
            }
          }
        }

        let grain_sum = 0.0;
        for (let ind_x = 0; ind_x < this.Nx; ind_x++) {
          for (let ind_y = 0; ind_y < this.Ny; ind_y++) {
            grain_sum += this.etas[ind_grain][ind_x][ind_y];
          }
        }

        grain_sum = grain_sum / (this.Nx * this.Ny);
        if (grain_sum < 0.001) {
          this.glist[ind_grain] = 0;
        }
      }
    }
  }

  convertToSaveData() {
    const data_out = Array(this.Nx)
      .fill()
      .map(() => Array(this.Ny).fill(0));

    for (let igrain = 0; igrain < this.ngrains; igrain++) {
      for (let i = 0; i < this.Nx; i++) {
        for (let j = 0; j < this.Ny; j++) {
          data_out[i][j] += this.etas[igrain][i][j] * this.etas[igrain][i][j];
        }
      }
    }
    return data_out;
  }

  runSimulation() {
    const grainCenters = this.generateRandomPoints(this.ngrains);
    this.generatePolycrystalline(grainCenters);
    const results = [];

    for (let ind_step = 0; ind_step < this.tEnd; ind_step++) {
      this.grainGrowthStep();

      if (ind_step % this.outInterval === 0) {
        let data_out = this.convertToSaveData();
        results.push({
          time: ind_step,
          data: data_out,
        });
        console.log(ind_step)
      }
    }
    return results;
  }
}
