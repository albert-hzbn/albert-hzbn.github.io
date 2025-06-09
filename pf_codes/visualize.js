const canvas = document.getElementById("gridCanvas");
const ctx = canvas.getContext("2d");
const slider = document.getElementById("timeSlider");
const timeStepDisplay = document.getElementById("timeStep");
const playButton = document.getElementById("playButton");
const resetButton = document.getElementById("resetButton");
const frameDelayInput = document.getElementById("frameDelay");
const startCalculationButton = document.getElementById("startCalculation");

// Define grid dimensions from canvas
const gridWidth = canvas.width;
const gridHeight = canvas.height;

// Animation configuration
let isPlaying = false;
let animationId = null;
let frameDelay = parseInt(frameDelayInput.value);
let lastFrameTime = 0;

function sleep(milliseconds) {
  const date = Date.now();
  let currentDate = null;
  do {
    currentDate = Date.now();
  } while (currentDate - date < milliseconds);
}

let simulationResults;

async function startCalculations(simulation_type, simulationParams) {
  console.log(simulationParams);
  console.log("Started");
  document.getElementById("status").textContent = "Running...";
  await new Promise((resolve) => setTimeout(resolve, 1000)); // Pause for 1 second
  console.log("now running");

  // Run simulation
  if (simulation_type == "spinodal") {
    const simulation = new CahnHilliardSimulation(simulationParams);
    simulationResults = simulation.runSimulation();
  } else if (simulation_type == "grain_growth") {
    const simulation = new GrainGrowthSimulation(simulationParams);
    simulationResults = simulation.runSimulation();
  }

  document.getElementById("status").textContent = "Completed";
  // console.log(simulationResults);

  // Set up slider and draw initial grid
  slider.max = simulationResults.length - 1;

  slider.addEventListener("input", (e) => {
    const timeStep = parseInt(e.target.value);
    timeStepDisplay.textContent = `Time Step: ${timeStep}`;
    drawGrid(simulationResults, timeStep);
  });

  // Initial draw
  drawGrid(simulationResults, 0);
}

function getColorForValue(value) {
  const intensity = Math.floor(value * 255);
  return `rgb(${intensity}, ${intensity}, ${intensity})`;
}

function drawGrid(data, timeStep) {
  ctx.clearRect(0, 0, gridWidth, gridHeight);

  const currentData = data[timeStep];
  const size_x = currentData.data.length;
  const size_y = currentData.data[0].length;

  const cellWidth = gridWidth / size_x;
  const cellHeight = gridHeight / size_y;

  const values2D = currentData.data;

  // Flatten the 2D array to compute min and max
  const flatValues = values2D.flat();  // or use .reduce if .flat isn't available
  const min = Math.min(...flatValues);
  const max = Math.max(...flatValues);

  // Normalize each value in the 2D array
  const normalizedData = (max === min)
    ? values2D.map(row => row.map(_ => 0.5))  // or 0 or 1, based on your preference
    : values2D.map(row => row.map(value => (value - min) / (max - min)));


  // Draw cells
  for (let x = 0; x < size_x; x++) {
    for (let y = 0; y < size_y; y++) {
      const value = normalizedData[x][y];
      ctx.fillStyle = getColorForValue(value);
      ctx.fillRect(x * cellWidth, y * cellHeight, cellWidth, cellHeight);
    }
  }
  console.log(currentData.data,normalizedData);

  // Grid lines
  ctx.strokeStyle = "#000000";
  ctx.lineWidth = 0.001;

  for (let x = 0; x <= size_x; x++) {
    ctx.beginPath();
    ctx.moveTo(x * cellWidth, 0);
    ctx.lineTo(x * cellWidth, gridHeight);
    ctx.stroke();
  }

  for (let y = 0; y <= size_y; y++) {
    ctx.beginPath();
    ctx.moveTo(0, y * cellHeight);
    ctx.lineTo(gridWidth, y * cellHeight);
    ctx.stroke();
  }
}

function animate(timestamp) {
  if (!isPlaying) return;

  if (timestamp - lastFrameTime >= frameDelay) {
    let currentStep = parseInt(slider.value);
    if (currentStep >= simulationResults.length - 1) {
      currentStep = 0;
    } else {
      currentStep++;
    }

    slider.value = currentStep;
    timeStepDisplay.textContent = `Time Step: ${currentStep}`;
    drawGrid(simulationResults, currentStep);

    lastFrameTime = timestamp;
  }

  animationId = requestAnimationFrame(animate);
}

function togglePlay() {
  isPlaying = !isPlaying;
  playButton.textContent = isPlaying ? "Pause" : "Play";

  if (isPlaying) {
    animationId = requestAnimationFrame(animate);
  } else {
    cancelAnimationFrame(animationId);
  }
}

function reset() {
  isPlaying = false;
  playButton.textContent = "Play";
  cancelAnimationFrame(animationId);
  slider.value = 0;
  timeStepDisplay.textContent = "Time Step: 0";
  drawGrid(simulationResults, 0);
}

frameDelayInput.addEventListener("change", (e) => {
  frameDelay = parseInt(e.target.value);
});

playButton.addEventListener("click", togglePlay);
resetButton.addEventListener("click", reset);
