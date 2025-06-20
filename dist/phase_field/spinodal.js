class CahnHilliardSimulation {
  constructor(params) {
    this.Nx = params.Nx;
    this.Ny = params.Ny;
    this.dx = params.dx;
    this.dy = params.dy;
    this.M = params.M;
    this.epsilon = params.epsilon;
    this.avgComp = params.avgComp;
    this.fluctuation = params.fluctuation;
    this.aVal = params.aVal;
    this.kappa = params.kappa;
    this.dt = params.dt;
    this.tEnd = params.tEnd;
    this.outInterval = params.outInterval;
  }

  zeros(rows, cols) {
    return Array(rows)
      .fill()
      .map(() => Array(cols).fill(0));
  }

  initialCondition() {
    const lowerBound = this.avgComp * (1 - this.fluctuation);
    const upperBound = this.avgComp * (1 + this.fluctuation);
    const comp = this.zeros(this.Nx, this.Ny);

    for (let indX = 0; indX < this.Nx; indX++) {
      for (let indY = 0; indY < this.Ny; indY++) {
        comp[indX][indY] =
          Math.random() * (upperBound - lowerBound) + lowerBound;
      }
    }
    return comp;
  }

  deepCopy(arr) {
    return arr.map((row) => [...row]);
  }

  applyPeriodicBC(array) {
    const result = this.deepCopy(array);
    const rows = array.length;
    const cols = array[0].length;

    for (let j = 0; j < cols; j++) {
      result[0][j] = array[rows - 2][j];
      result[rows - 1][j] = array[1][j];
    }

    for (let i = 0; i < rows; i++) {
      result[i][0] = array[i][cols - 2];
      result[i][cols - 1] = array[i][1];
    }

    return result;
  }

  cahnHilliardStep(comp) {
    const mu = this.zeros(this.Nx, this.Ny);
    const compExt = this.applyPeriodicBC(comp);

    for (let indX = 1; indX < this.Nx - 1; indX++) {
      for (let indY = 1; indY < this.Ny - 1; indY++) {
        const laplaceComp =
          (compExt[indX - 1][indY] -
            2 * compExt[indX][indY] +
            compExt[indX + 1][indY]) /
            (this.dx * this.dx) +
          (compExt[indX][indY - 1] -
            2 * compExt[indX][indY] +
            compExt[indX][indY + 1]) /
            (this.dy * this.dy);

        const c = compExt[indX][indY];
        const dfDc = 2 * this.aVal * c * (1 - c) * (1 - 2 * c);

        mu[indX][indY] = dfDc - this.kappa * laplaceComp;
      }
    }

    const muExt = this.applyPeriodicBC(mu);

    for (let indX = 1; indX < this.Nx - 1; indX++) {
      for (let indY = 1; indY < this.Ny - 1; indY++) {
        const laplaceMu =
          (muExt[indX - 1][indY] -
            2 * muExt[indX][indY] +
            muExt[indX + 1][indY]) /
            (this.dx * this.dx) +
          (muExt[indX][indY - 1] -
            2 * muExt[indX][indY] +
            muExt[indX][indY + 1]) /
            (this.dy * this.dy);

        compExt[indX][indY] = comp[indX][indY] + this.dt * this.M * laplaceMu;
      }
    }

    return compExt;
  }

  runSimulation() {
    let comp = this.initialCondition();
    const results = [];

    for (let ind_step = 0; ind_step < this.tEnd; ind_step++) {
      comp = this.cahnHilliardStep(comp);
      if (ind_step % this.outInterval === 0) {
        results.push({
          time: ind_step,
          data: this.deepCopy(comp),
        });
      }
      console.log(ind_step)
    }

    return results;
  }
}
