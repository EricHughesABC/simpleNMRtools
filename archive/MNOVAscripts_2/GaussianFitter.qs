



/**
 * Nelder-Mead Simplex Algorithm for Gaussian Peak Fitting
 * ECMAScript 5 (ECMA-262) Compatible
 */

function GaussianFitter() {
    'use strict';
    
    // Gaussian function: f(x) = a * exp(-((x - b)^2) / (2 * c^2)) + d
    // where: a = amplitude, b = center, c = sigma (width), d = baseline
    function gaussian(x, params) {
        var a = params[0]; // amplitude
        var b = params[1]; // center
        var c = params[2]; // sigma
        var d = params[3]; // baseline offset
        return a * Math.exp(-Math.pow(x - b, 2) / (2 * c * c)) + d;
    }
    
    // Calculate sum of squared residuals
    function calculateError(params, xData, yData) {
        var error = 0;
        var i, predicted, residual;
        
        for (i = 0; i < xData.length; i++) {
            predicted = gaussian(xData[i], params);
            residual = yData[i] - predicted;
            error += residual * residual;
        }
        
        return error;
    }
    
    // Nelder-Mead Simplex Algorithm
    function nelderMead(objective, initialGuess, options) {
        options = options || {};
        var maxIterations = options.maxIterations || 1000;
        var tolerance = options.tolerance || 1e-8;
        var alpha = options.alpha || 1.0;  // reflection
        var gamma = options.gamma || 2.0;  // expansion
        var rho = options.rho || 0.5;      // contraction
        var sigma = options.sigma || 0.5;  // shrink
        
        var n = initialGuess.length;
        var simplex = [];
        var i, j;
        
        // Initialize simplex
        simplex[0] = initialGuess.slice();
        for (i = 1; i <= n; i++) {
            simplex[i] = initialGuess.slice();
            simplex[i][i - 1] += simplex[i][i - 1] * 0.05 + 0.00025;
        }
        
        // Evaluate initial simplex
        var values = [];
        for (i = 0; i <= n; i++) {
            values[i] = objective(simplex[i]);
        }
        
        var iteration = 0;
        var centroid, reflected, expanded, contracted, shrunk;
        var reflectedValue, expandedValue, contractedValue;
        var worstIndex, secondWorstIndex, bestIndex;
        
        while (iteration < maxIterations) {
            // Sort simplex by function values
            var indices = [];
            for (i = 0; i <= n; i++) {
                indices[i] = i;
            }
            
            indices.sort(function(a, b) {
                return values[a] - values[b];
            });
            
            bestIndex = indices[0];
            secondWorstIndex = indices[n - 1];
            worstIndex = indices[n];
            
            // Check convergence
            var range = values[worstIndex] - values[bestIndex];
            if (range < tolerance) {
                break;
            }
            
            // Calculate centroid (excluding worst point)
            centroid = new Array(n);
            for (i = 0; i < n; i++) {
                centroid[i] = 0;
                for (j = 0; j < n; j++) {
                    centroid[i] += simplex[indices[j]][i];
                }
                centroid[i] /= n;
            }
            
            // Reflection
            reflected = new Array(n);
            for (i = 0; i < n; i++) {
                reflected[i] = centroid[i] + alpha * (centroid[i] - simplex[worstIndex][i]);
            }
            reflectedValue = objective(reflected);
            
            if (reflectedValue >= values[bestIndex] && reflectedValue < values[secondWorstIndex]) {
                // Accept reflection
                simplex[worstIndex] = reflected;
                values[worstIndex] = reflectedValue;
            } else if (reflectedValue < values[bestIndex]) {
                // Try expansion
                expanded = new Array(n);
                for (i = 0; i < n; i++) {
                    expanded[i] = centroid[i] + gamma * (reflected[i] - centroid[i]);
                }
                expandedValue = objective(expanded);
                
                if (expandedValue < reflectedValue) {
                    simplex[worstIndex] = expanded;
                    values[worstIndex] = expandedValue;
                } else {
                    simplex[worstIndex] = reflected;
                    values[worstIndex] = reflectedValue;
                }
            } else {
                // Try contraction
                contracted = new Array(n);
                for (i = 0; i < n; i++) {
                    contracted[i] = centroid[i] + rho * (simplex[worstIndex][i] - centroid[i]);
                }
                contractedValue = objective(contracted);
                
                if (contractedValue < values[worstIndex]) {
                    simplex[worstIndex] = contracted;
                    values[worstIndex] = contractedValue;
                } else {
                    // Shrink simplex
                    for (i = 0; i <= n; i++) {
                        if (i !== bestIndex) {
                            for (j = 0; j < n; j++) {
                                simplex[i][j] = simplex[bestIndex][j] + 
                                    sigma * (simplex[i][j] - simplex[bestIndex][j]);
                            }
                            values[i] = objective(simplex[i]);
                        }
                    }
                }
            }
            
            iteration++;
        }
        
        // Find best solution
        var minValue = values[0];
        var minIndex = 0;
        for (i = 1; i <= n; i++) {
            if (values[i] < minValue) {
                minValue = values[i];
                minIndex = i;
            }
        }
        
        return {
            solution: simplex[minIndex],
            value: minValue,
            iterations: iteration,
            converged: iteration < maxIterations
        };
    }
    
    // Generate initial guess from data
    function generateInitialGuess(xData, yData) {
        var maxY = Math.max.apply(Math, yData);
        var minY = Math.min.apply(Math, yData);
        var maxIndex = yData.indexOf(maxY);
        var centerX = xData[maxIndex]-10; // Center at peak position with some offset
        
        // Estimate baseline as minimum value or average of first/last few points
        var n = Math.min(5, Math.floor(yData.length / 4)); // Use 5 points or 25% of data, whichever is smaller
        var baselineSum = 0;
        var i;
        
        for (i = 0; i < n; i++) {
            baselineSum += yData[i] + yData[yData.length - 1 - i];
        }
        var estimatedBaseline = baselineSum / (2 * n);
        
        // Amplitude is peak height above baseline
        var estimatedAmplitude = maxY - estimatedBaseline;
        
        // Estimate width as 1/4 of data range
        var minX = Math.min.apply(Math, xData);
        var maxX = Math.max.apply(Math, xData);
        // var estimatedWidth = (maxX - minX) / 10;
        var estimatedWidth = 0.9;

        print('\nGenerating initial guess for Gaussian parameters:\n');

        print('Initial guess: Amplitude = ' + estimatedAmplitude.toFixed(2));
        print('Center = ' + centerX.toFixed(2));
        print('Width (Sigma) = ' + estimatedWidth.toFixed(2));
        print('Baseline = ' + estimatedBaseline.toFixed(2) + '\n');
        
        return [estimatedAmplitude, centerX, estimatedWidth, estimatedBaseline];
    }
    
    // Generate fitted curve coordinates
    function generateFittedCurve(fittedParams, xData, options) {
        options = options || {};
        var numPoints = options.numPoints || (xData.length * 2); // Higher resolution by default
        var xMin = options.xMin || Math.min.apply(Math, xData);
        var xMax = options.xMax || Math.max.apply(Math, xData);
        
        var xFitted = [];
        var yFitted = [];
        var i, x;
        
        for (i = 0; i < numPoints; i++) {
            x = xMin + (xMax - xMin) * i / (numPoints - 1);
            xFitted.push(x);
            yFitted.push(gaussian(x, fittedParams));
        }
        
        return {
            x: xFitted,
            y: yFitted
        };
    }
    
    // Main fitting function
    function fitGaussian(xData, yData, options) {
        options = options || {};
        var initialGuess = options.initialGuess || generateInitialGuess(xData, yData);
        
        // Create objective function
        function objective(params) {
            // Add constraints to prevent invalid parameters
            if (params[0] <= 0 || params[2] <= 0) {
                return 1e10; // Large penalty for invalid parameters (amplitude and sigma must be positive)
            }
            return calculateError(params, xData, yData);
        }
        
        var result = nelderMead(objective, initialGuess, options);
        
        // Calculate R-squared
        var totalSumSquares = 0;
        var residualSumSquares = result.value;
        var meanY = 0;
        var i;
        
        for (i = 0; i < yData.length; i++) {
            meanY += yData[i];
        }
        meanY /= yData.length;
        
        for (i = 0; i < yData.length; i++) {
            totalSumSquares += Math.pow(yData[i] - meanY, 2);
        }
        
        var rSquared = 1 - (residualSumSquares / totalSumSquares);
        
        // Generate fitted curve at original x points
        var fittedAtOriginal = [];
        for (i = 0; i < xData.length; i++) {
            fittedAtOriginal.push(gaussian(xData[i], result.solution));
        }
        
        return {
            amplitude: result.solution[0],
            center: result.solution[1],
            sigma: result.solution[2],
            baseline: result.solution[3],
            rSquared: rSquared,
            error: result.value,
            iterations: result.iterations,
            converged: result.converged,
            fittedY: fittedAtOriginal, // Y values at original X points
            gaussianFunction: function(x) {
                return gaussian(x, result.solution);
            },
            generateCurve: function(options) {
                return generateFittedCurve(result.solution, xData, options);
            }
        };
    }
    
    // Public API
    return {
        fit: fitGaussian,
        gaussian: gaussian,
        generateInitialGuess: generateInitialGuess,
        // Optional: fit without baseline (3-parameter version)
        fitNoBaseline: function(xData, yData, options) {
            // Wrapper that fixes baseline to 0
            var oldGaussian = gaussian;
            gaussian = function(x, params) {
                return oldGaussian(x, [params[0], params[1], params[2], 0]);
            };
            
            var result = fitGaussian(xData, yData, options);
            gaussian = oldGaussian; // Restore original function
            
            // Remove baseline from result for clarity
            delete result.baseline;
            return result;
        }
    };
}

// Test Data Generator
function TestDataGenerator() {
    'use strict';
    
    // Generate random number with normal distribution (Box-Muller transform)
    function normalRandom(mean, stdDev) {
        mean = mean || 0;
        stdDev = stdDev || 1;
        
        var u1 = Math.random();
        var u2 = Math.random();
        var z0 = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
        return z0 * stdDev + mean;
    }
    
    // Generate basic Gaussian peak with noise
    function generateGaussianPeak(params) {
        params = params || {};
        var amplitude = params.amplitude || 5.0;
        var center = params.center || 0.0;
        var sigma = params.sigma || 1.0;
        var numPoints = params.numPoints || 100;
        var xMin = params.xMin || (center - 4 * sigma);
        var xMax = params.xMax || (center + 4 * sigma);
        var noiseLevel = params.noiseLevel || 0.1;
        var baseline = params.baseline || 0.0;
        
        var xData = [];
        var yData = [];
        var i, x, signal, noise;
        
        for (i = 0; i < numPoints; i++) {
            x = xMin + (xMax - xMin) * i / (numPoints - 1);
            signal = amplitude * Math.exp(-Math.pow(x - center, 2) / (2 * sigma * sigma));
            noise = normalRandom(0, noiseLevel * amplitude);
            
            xData.push(x);
            yData.push(signal + baseline + noise);
        }
        
        return {
            x: xData,
            y: yData,
            trueParameters: {
                amplitude: amplitude,
                center: center,
                sigma: sigma,
                baseline: baseline,
                noiseLevel: noiseLevel
            }
        };
    }
    
    // Generate multiple overlapping peaks
    function generateMultiplePeaks(peaks) {
        peaks = peaks || [
            { amplitude: 3, center: -2, sigma: 0.8 },
            { amplitude: 5, center: 1, sigma: 1.2 },
            { amplitude: 2, center: 4, sigma: 0.6 }
        ];
        
        var numPoints = 150;
        var xMin = -6;
        var xMax = 7;
        var noiseLevel = 0.05;
        var baseline = 0.2;
        
        var xData = [];
        var yData = [];
        var i, j, x, totalSignal, noise;
        
        for (i = 0; i < numPoints; i++) {
            x = xMin + (xMax - xMin) * i / (numPoints - 1);
            totalSignal = baseline;
            
            for (j = 0; j < peaks.length; j++) {
                totalSignal += peaks[j].amplitude * 
                    Math.exp(-Math.pow(x - peaks[j].center, 2) / (2 * peaks[j].sigma * peaks[j].sigma));
            }
            
            noise = normalRandom(0, noiseLevel * Math.max.apply(Math, peaks.map(function(p) { return p.amplitude; })));
            
            xData.push(x);
            yData.push(totalSignal + noise);
        }
        
        return {
            x: xData,
            y: yData,
            truePeaks: peaks,
            baseline: baseline,
            noiseLevel: noiseLevel
        };
    }
    
    // Generate peak with linear background
    function generatePeakWithBackground(params) {
        params = params || {};
        var amplitude = params.amplitude || 4.0;
        var center = params.center || 2.0;
        var sigma = params.sigma || 1.5;
        var slope = params.slope || 0.3;
        var intercept = params.intercept || 1.0;
        var numPoints = params.numPoints || 80;
        var xMin = params.xMin || -3;
        var xMax = params.xMax || 7;
        var noiseLevel = params.noiseLevel || 0.08;
        
        var xData = [];
        var yData = [];
        var i, x, signal, background, noise;
        
        for (i = 0; i < numPoints; i++) {
            x = xMin + (xMax - xMin) * i / (numPoints - 1);
            signal = amplitude * Math.exp(-Math.pow(x - center, 2) / (2 * sigma * sigma));
            background = slope * x + intercept;
            noise = normalRandom(0, noiseLevel * amplitude);
            
            xData.push(x);
            yData.push(signal + background + noise);
        }
        
        return {
            x: xData,
            y: yData,
            trueParameters: {
                amplitude: amplitude,
                center: center,
                sigma: sigma,
                slope: slope,
                intercept: intercept,
                noiseLevel: noiseLevel
            }
        };
    }
    
    // Generate asymmetric peak (skewed Gaussian)
    function generateAsymmetricPeak(params) {
        params = params || {};
        var amplitude = params.amplitude || 6.0;
        var center = params.center || 1.0;
        var sigma = params.sigma || 1.0;
        var skewness = params.skewness || 2.0;
        var numPoints = params.numPoints || 120;
        var xMin = params.xMin || -4;
        var xMax = params.xMax || 6;
        var noiseLevel = params.noiseLevel || 0.12;
        
        var xData = [];
        var yData = [];
        var i, x, t, signal, noise;
        
        for (i = 0; i < numPoints; i++) {
            x = xMin + (xMax - xMin) * i / (numPoints - 1);
            t = (x - center) / sigma;
            
            // Skewed normal approximation
            signal = amplitude * Math.exp(-t * t / 2) * (1 + Math.sign(t) * Math.tanh(skewness * Math.abs(t) / 2));
            noise = normalRandom(0, noiseLevel * amplitude);
            
            xData.push(x);
            yData.push(Math.max(0, signal + noise)); // Ensure non-negative
        }
        
        return {
            x: xData,
            y: yData,
            trueParameters: {
                amplitude: amplitude,
                center: center,
                sigma: sigma,
                skewness: skewness,
                noiseLevel: noiseLevel
            }
        };
    }
    
    // Generate test datasets collection
    function generateTestSuite() {
        return {
            // Clean, well-behaved peak
            clean: generateGaussianPeak({
                amplitude: 10,
                center: 0,
                sigma: 1.5,
                numPoints: 50,
                noiseLevel: 0.02
            }),
            
            // Noisy peak
            noisy: generateGaussianPeak({
                amplitude: 3,
                center: 2,
                sigma: 0.8,
                numPoints: 60,
                noiseLevel: 0.25
            }),
            
            // Narrow peak
            narrow: generateGaussianPeak({
                amplitude: 8,
                center: -1,
                sigma: 0.3,
                numPoints: 80,
                noiseLevel: 0.1
            }),
            
            // Broad peak
            broad: generateGaussianPeak({
                amplitude: 2,
                center: 1,
                sigma: 3,
                numPoints: 100,
                noiseLevel: 0.05
            }),
            
            // Peak with baseline
            withBaseline: generateGaussianPeak({
                amplitude: 4,
                center: 0.5,
                sigma: 1.2,
                baseline: 2.5,
                numPoints: 70,
                noiseLevel: 0.15
            }),
            
            // Multiple peaks (challenging for single peak fitting)
            multiplePeaks: generateMultiplePeaks(),
            
            // Peak with linear background
            withBackground: generatePeakWithBackground(),
            
            // Asymmetric peak
            asymmetric: generateAsymmetricPeak()
        };
    }
    
    return {
        generateGaussianPeak: generateGaussianPeak,
        generateMultiplePeaks: generateMultiplePeaks,
        generatePeakWithBackground: generatePeakWithBackground,
        generateAsymmetricPeak: generateAsymmetricPeak,
        generateTestSuite: generateTestSuite
    };
}

// Example usage and testing
function runTests() {
    'use strict';
    
    var generator = TestDataGenerator();
    var fitter = GaussianFitter();
    var testSuite = generator.generateTestSuite();
    
    print('=== Gaussian Peak Fitting Test Results ===\n');
    
    // Test each dataset
    var testNames = ['clean', 'noisy', 'narrow', 'broad', 'withBaseline'];
    var i, testName, data, result, trueParams;
    
    for (i = 0; i < testNames.length; i++) {
        testName = testNames[i];
        data = testSuite[testName];
        result = fitter.fit(data.x, data.y);
        trueParams = data.trueParameters;
        
        print('Test: ' + testName.toUpperCase());
        print('True params  - Amp: ' + trueParams.amplitude.toFixed(2) + 
                   ', Center: ' + trueParams.center.toFixed(2) + 
                   ', Sigma: ' + trueParams.sigma.toFixed(2) +
                   ', Baseline: ' + (trueParams.baseline || 0).toFixed(2));
        print('Fitted params - Amp: ' + result.amplitude.toFixed(2) + 
                   ', Center: ' + result.center.toFixed(2) + 
                   ', Sigma: ' + result.sigma.toFixed(2) +
                   ', Baseline: ' + result.baseline.toFixed(2));
        print('R-squared: ' + result.rSquared.toFixed(4) + 
                   ', Iterations: ' + result.iterations + 
                   ', Converged: ' + result.converged);
        print('---');
    }
    
    return testSuite;
}

// Generate specific test case
function generateSpecificTest() {
    var generator = TestDataGenerator();
    
    // Create a clean test case
    return generator.generateGaussianPeak({
        amplitude: 7.5,
        center: 1.2,
        sigma: 0.8,
        numPoints: 200,
        xMin: -5,
        xMax: 7,
        noiseLevel: 0.01,
        baseline: 2.0
    });
}

// Plotting and visualization helpers
function PlottingHelpers() {
    'use strict';
    
    // Generate CSV output for external plotting
    function generateCSV(originalData, fittedResult, options) {
        options = options || {};
        var highResolution = options.highResolution !== false; // Default true
        
        var lines = ['X,Y_Original,Y_Fitted'];
        var i;
        
        if (highResolution) {
            // Generate high-resolution fitted curve
            var curve = fittedResult.generateCurve({ numPoints: 200 });
            
            // Add original data points
            for (i = 0; i < originalData.x.length; i++) {
                lines.push(originalData.x[i] + ',' + originalData.y[i] + ',' + fittedResult.fittedY[i]);
            }
            
            // Add separator and high-res curve
            lines.push('# High-resolution fitted curve');
            for (i = 0; i < curve.x.length; i++) {
                lines.push(curve.x[i] + ',,' + curve.y[i]);
            }
        } else {
            // Just original points with fitted values
            for (i = 0; i < originalData.x.length; i++) {
                lines.push(originalData.x[i] + ',' + originalData.y[i] + ',' + fittedResult.fittedY[i]);
            }
        }
        
        return lines.join('\n');
    }
    
    // Generate simple text table
    function generateTable(originalData, fittedResult) {
        var lines = [];
        var i, x, yOrig, yFit, residual;
        
        lines.push('    X        Y_Orig    Y_Fitted   Residual');
        lines.push('----------------------------------------');
        
        for (i = 0; i < originalData.x.length; i++) {
            x = originalData.x[i];
            yOrig = originalData.y[i];
            yFit = fittedResult.fittedY[i];
            residual = yOrig - yFit;
            
            lines.push(
                x.toFixed(3).padStart(8) + '  ' +
                yOrig.toFixed(3).padStart(8) + '  ' +
                yFit.toFixed(3).padStart(8) + '  ' +
                residual.toFixed(3).padStart(8)
            );
        }
        
        return lines.join('\n');
    }
    
    // Generate plotting instructions for common tools
    function generatePlottingInstructions(dataVarName) {
        dataVarName = dataVarName || 'result';
        
        return {
            matplotlib: [
                '# Python/Matplotlib plotting:',
                'import matplotlib.pyplot as plt',
                'import numpy as np',
                '',
                '# Original data',
                'x_orig = ' + dataVarName + '.original.x',
                'y_orig = ' + dataVarName + '.original.y',
                '',
                '# High-resolution fitted curve',
                'curve = ' + dataVarName + '.fitted.generateCurve({numPoints: 200})',
                'x_fit = curve.x',
                'y_fit = curve.y',
                '',
                'plt.figure(figsize=(10, 6))',
                'plt.scatter(x_orig, y_orig, color="blue", alpha=0.7, label="Original Data")',
                'plt.plot(x_fit, y_fit, color="red", linewidth=2, label="Fitted Gaussian")',
                'plt.xlabel("X")',
                'plt.ylabel("Y")',
                'plt.legend()',
                'plt.grid(True, alpha=0.3)',
                'plt.title("Gaussian Peak Fit")',
                'plt.show()'
            ].join('\n'),
            
            excel: [
                '# Excel plotting:',
                '1. Copy the CSV data generated by generateCSV()',
                '2. Paste into Excel',
                '3. Select all data',
                '4. Insert > Scatter Chart',
                '5. Format series to show points for original data and smooth line for fitted curve'
            ].join('\n'),
            
            gnuplot: [
                '# Gnuplot commands:',
                'set xlabel "X"',
                'set ylabel "Y"',
                'set grid',
                'plot "data.csv" using 1:2 with points title "Original Data", \\',
                '     "data.csv" using 1:3 with lines title "Fitted Gaussian"'
            ].join('\n')
        };
    }
    
    return {
        generateCSV: generateCSV,
        generateTable: generateTable,
        generatePlottingInstructions: generatePlottingInstructions
    };
}

// Complete example with plotting data
function completeExample() {
    'use strict';
    
    var generator = TestDataGenerator();
    var fitter = GaussianFitter();
    var plotter = PlottingHelpers();
    
    // Generate test data
    var testData = generator.generateGaussianPeak({
        amplitude: 5.0,
        center: 1.5,
        sigma: 1.2,
        numPoints: 30,
        noiseLevel: 0.15,
        baseline: 0.5
    });
    
    // Fit Gaussian
    var result = fitter.fit(testData.x, testData.y);
    
    // Print results
    print('=== Gaussian Fitting Results ===');
    print('True parameters:');
    print('  Amplitude: ' + testData.trueParameters.amplitude);
    print('  Center: ' + testData.trueParameters.center);
    print('  Sigma: ' + testData.trueParameters.sigma);
    print('');
    print('Fitted parameters:');
    print('  Amplitude: ' + result.amplitude.toFixed(4));
    print('  Center: ' + result.center.toFixed(4));
    print('  Sigma: ' + result.sigma.toFixed(4));
    print('  Baseline: ' + result.baseline.toFixed(4));
    print('  R-squared: ' + result.rSquared.toFixed(4));
    print('');
    
    // Method 1: Get Y coordinates at original X points
    print('Method 1 - Y values at original X points:');
    print('result.fittedY contains:', result.fittedY.length, 'points');
    
    // Method 2: Evaluate function at any X value
    print('Method 2 - Evaluate at specific point:');
    print('f(2.0) =', result.gaussianFunction(2.0).toFixed(4));
    
    // Method 3: Generate high-resolution curve
    var highResCurve = result.generateCurve({ numPoints: 100 });
    print('Method 3 - High-resolution curve:');
    print('Generated', highResCurve.x.length, 'points for smooth plotting');
    
    // Method 4: Generate CSV for external plotting
    var csvData = plotter.generateCSV(testData, result);
    print('Method 4 - CSV data generated (', csvData.split('\n').length, 'lines)');
    
    // Method 5: Generate comparison table
    print('Method 5 - Data comparison table:');
    print(plotter.generateTable(testData, result));
    
    return {
        original: testData,
        fitted: result,
        highResCurve: highResCurve,
        csvData: csvData,
        plottingInstructions: plotter.generatePlottingInstructions('myResult')
    };
}

// Usage examples for getting fitted curve coordinates:
function usageExamples() {
    'use strict';
    
    print('=== Usage Examples for Plotting ===\n');
    
    // Example setup
    var xData = [-2, -1, 0, 1, 2, 3, 4];
    var yData = [0.5, 1.2, 3.8, 5.1, 3.9, 1.1, 0.4];
    var fitter = GaussianFitter();
    var result = fitter.fit(xData, yData);
    
    print('1. Get fitted Y values at original X points:');
    print('   var fittedY = result.fittedY;');
    print('   // fittedY.length === xData.length');
    print('   Values:', result.fittedY.map(function(y) { return y.toFixed(2); }));
    
    print('\n2. Evaluate fitted function at any X value:');
    print('   var y_at_2_5 = result.gaussianFunction(2.5);');
    print('   Result:', result.gaussianFunction(2.5).toFixed(4));
    
    print('\n3. Generate high-resolution curve for smooth plotting:');
    print('   var curve = result.generateCurve({ numPoints: 50 });');
    print('   // curve.x and curve.y arrays with 50 points each');
    var curve = result.generateCurve({ numPoints: 10 }); // Smaller for display
    print('   Sample X:', curve.x.slice(0, 5).map(function(x) { return x.toFixed(2); }));
    print('   Sample Y:', curve.y.slice(0, 5).map(function(y) { return y.toFixed(2); }));
    
    print('\n4. Generate curve over custom X range:');
    print('   var customCurve = result.generateCurve({');
    print('       numPoints: 100,');
    print('       xMin: -3,');
    print('       xMax: 5');
    print('   });');
    
    print('\n5. For loop to generate your own points:');
    print('   for (var i = 0; i < numPoints; i++) {');
    print('       var x = xMin + (xMax - xMin) * i / (numPoints - 1);');
    print('       var y = result.gaussianFunction(x);');
    print('       // Use x, y for plotting');
    print('   }');
}

// Uncomment to run examples:
// var example = completeExample();
// usageExamples();


function simplex_eeh() {
    'use strict';


var yData = [-66.177734375,113.07421875,-2.564453125,-35.474609375,11.798828125,-136.345703125,-178.1953125,-160.4765625,-367.9140625,-853.85546875,14576.62890625,1248.9921875,14550.0078125,2770.828125,602.162109375,373.736328125,228.67578125,21.322265625,135.46484375,26.25]
var xData = [141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160]
    // Initialize Gaussian fitter

   
    var fitter = GaussianFitter();
    var initialGuess = fitter.generateInitialGuess(xData, yData);

     print('Initial Guess:\n'+ initialGuess);
 

    var fitResult = fitter.fit(xData, yData);

    var initialPlot = fitResult.generateCurve(initialGuess, xData);
    
    // // Generate specific test data
    // var testData = generateSpecificTest();
    
    // // Fit the Gaussian peak
    // var fitResult = fitter.fit(testData.x, testData.y);
    
    // Print results
    print('Fitted Gaussian Peak:');
    print('Amplitude: ' + fitResult.amplitude.toFixed(2));
    print('Center: ' + fitResult.center.toFixed(2));
    print('Sigma: ' + fitResult.sigma.toFixed(2));
    print('Baseline: ' + fitResult.baseline.toFixed(2));
    print('R-squared: ' + fitResult.rSquared.toFixed(4));
    print('Iterations: ' + fitResult.iterations);
    print('Converged: ' + fitResult.converged);
    print('xxx = ' + JSON.stringify(xData));
    print('yyy = ' + JSON.stringify(yData));
    print('yyy_res = ' + JSON.stringify(fitResult.fittedY));
    print('init_plot_x = ' + JSON.stringify(initialPlot.x));
    print('init_plot_y = ' + JSON.stringify(initialPlot.y));

    // print(fitResult);
}
// // Uncomment to run tests:
// // var testResults = runTests();
// // var specificTest = generateSpecificTest();