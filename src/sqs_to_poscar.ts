// Function for matrix multiplication
function matMul(matrix1: number[][], matrix2: number[][]): number[][] {
    const rows1 = matrix1.length;
    const cols1 = matrix1[0].length;
    const rows2 = matrix2.length;
    const cols2 = matrix2[0].length;

    if (cols1 !== rows2) {
        throw new Error('Cannot multiply matrices: dimensions do not match');
    }

    const result: number[][] = [];
    for (let i = 0; i < rows1; i++) {
        result[i] = [];
        for (let j = 0; j < cols2; j++) {
            let sum = 0;
            for (let k = 0; k < cols1; k++) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}


function convert_sqs_to_poscar() {
    const bestsqs = (<HTMLInputElement>document.getElementById("bestsqs")).value.trim();
    const bestsqsArr = bestsqs.split("\n");
    const bestsqsLen = bestsqsArr.length;
    const atomicPositions = bestsqsArr.slice(6, bestsqsLen);

    // Gets the elements prensent in the alloy 
    let elemArr: string[] = []
    for (let index = 0; index < atomicPositions.length; index++) {
        const atomPosStrArr = atomicPositions[index].split(" ");
        const elem: string = atomPosStrArr.slice(3)[0];
        if (!elemArr.includes(elem)) {
            elemArr.push(elem)
        }
    }

    // Get the atomic coordinates arragned and number of atoms for each element
    let countTotAtoms = 0;
    let numAtomsArr: number[] = [];
    let elemCoordArr: number[][] = [[]];
    for (let elemIndex = 0; elemIndex < elemArr.length; elemIndex++) {
        let countAtoms = 0;
        for (let index = 0; index < atomicPositions.length; index++) {
            if (elemArr[elemIndex] === atomicPositions[index].split(" ")[3]) {
                const atomPosStrArr = atomicPositions[index].split(" ");
                const elemCoordStr: string[] = atomPosStrArr.slice(0, 3);
                elemCoordArr[countTotAtoms].push(parseFloat(elemCoordStr[0]))
                elemCoordArr[countTotAtoms].push(parseFloat(elemCoordStr[1]))
                elemCoordArr[countTotAtoms].push(parseFloat(elemCoordStr[2]))
                elemCoordArr.push([]);
                elemCoordArr[countAtoms].push();
                countAtoms++;
                countTotAtoms++;
            }
        }
        numAtomsArr.push(countAtoms)
    }

    // Extracts unit cell and lattice vectors
    const unitCellStr = bestsqsArr.slice(0, 3);
    const latticeVectorsStr = bestsqsArr.slice(3, 6);
    let unitCellArr: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    let latticeVectorArr: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let index = 0; index < 3; index++) {
        const c1 = unitCellStr[index].split(" ")
        const c2 = latticeVectorsStr[index].split(" ")
        for (let index2 = 0; index2 < 3; index2++) {
            unitCellArr[index][index2] = parseFloat(c1[index2]);
            latticeVectorArr[index][index2] = parseFloat(c2[index2]);
        }
    }

    // Perform the matrix multiplication to obtain lattice vectors and atomic positions
    const newLatticeVectors = matMul(latticeVectorArr, unitCellArr);
    const newAtomPos = matMul(elemCoordArr, unitCellArr);

    // Converts data to POSCAR format
    let strPOSCAR = "POSCAR\n1.0\n";
    for (let index = 0; index < 3; index++) {
        strPOSCAR += `  ${newLatticeVectors[index][0].toFixed(10)}  ${newLatticeVectors[index][1].toFixed(10)}  ${newLatticeVectors[index][2].toFixed(10)}\n`;
    }
    strPOSCAR += "  ";
    for (let index = 0; index < elemArr.length; index++) {
        strPOSCAR += `${elemArr[index]} `;
    }
    strPOSCAR += "\n  ";
    for (let index = 0; index < elemArr.length; index++) {
        strPOSCAR += `${numAtomsArr[index]} `;
    }
    strPOSCAR += "\nCartesian\n";
    for (let index = 0; index < atomicPositions.length; index++) {
        strPOSCAR += `  ${newAtomPos[index][0].toFixed(10)}  ${newAtomPos[index][1].toFixed(10)}  ${newAtomPos[index][2].toFixed(10)}\n`;
    }

    // Writes POSCAR string to textarea
    (<HTMLInputElement>document.getElementById("poscar")).value = strPOSCAR;
}