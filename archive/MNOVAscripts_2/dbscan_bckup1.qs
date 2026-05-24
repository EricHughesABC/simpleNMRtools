

// // // DBSCAN Clustering Implementation - ES5 Compatible
// // function DBSCAN(eps, minPts) {
// //   this.eps = eps || 0.5;        // Maximum distance between two points to be neighbors
// //   this.minPts = minPts || 5;    // Minimum number of points to form a cluster
// // }

// // // Calculate Euclidean distance between two points
// // DBSCAN.prototype.distance = function(point1, point2) {
// //   var sum = 0;
// //   for (var i = 0; i < point1.length; i++) {
// //     sum += Math.pow(point1[i] - point2[i], 2);
// //   }
// //   return Math.sqrt(sum);
// // };

// // // Find all neighbors within eps distance
// // DBSCAN.prototype.regionQuery = function(dataset, pointIndex) {
// //   var neighbors = [];
// //   var point = dataset[pointIndex];
  
// //   for (var i = 0; i < dataset.length; i++) {
// //     if (this.distance(point, dataset[i]) <= this.eps) {
// //       neighbors.push(i);
// //     }
// //   }
// //   return neighbors;
// // };

// // // Main DBSCAN clustering algorithm
// // DBSCAN.prototype.fit = function(dataset) {
// //   var n = dataset.length;
// //   var labels = [];
// //   var clusterId = 0;
// //   var i, j;

// //   // Initialize labels array with -1 (unclassified)
// //   for (i = 0; i < n; i++) {
// //     labels[i] = -1;
// //   }

// //   for (i = 0; i < n; i++) {
// //     // Skip if point is already processed
// //     if (labels[i] !== -1) continue;

// //     // Find neighbors
// //     var neighbors = this.regionQuery(dataset, i);

// //     // Mark as noise if not enough neighbors
// //     if (neighbors.length < this.minPts) {
// //       labels[i] = -2; // Noise point
// //       continue;
// //     }

// //     // Start new cluster
// //     clusterId++;
// //     labels[i] = clusterId;

// //     // Process all neighbors
// //     var seedSet = neighbors.slice(); // Copy array
// //     j = 0;

// //     while (j < seedSet.length) {
// //       var currentPoint = seedSet[j];

// //       // Change noise to border point
// //       if (labels[currentPoint] === -2) {
// //         labels[currentPoint] = clusterId;
// //       }

// //       // Skip if already processed
// //       if (labels[currentPoint] !== -1) {
// //         j++;
// //         continue;
// //       }

// //       // Mark as part of cluster
// //       labels[currentPoint] = clusterId;

// //       // Find neighbors of current point
// //       var currentNeighbors = this.regionQuery(dataset, currentPoint);

// //       // If enough neighbors, add them to seed set
// //       if (currentNeighbors.length >= this.minPts) {
// //         for (var k = 0; k < currentNeighbors.length; k++) {
// //           var neighbor = currentNeighbors[k];
// //           var found = false;
          
// //           // Check if neighbor already in seedSet
// //           for (var m = 0; m < seedSet.length; m++) {
// //             if (seedSet[m] === neighbor) {
// //               found = true;
// //               break;
// //             }
// //           }
          
// //           if (!found) {
// //             seedSet.push(neighbor);
// //           }
// //         }
// //       }

// //       j++;
// //     }
// //   }

// //   return labels;
// // };

// // // Calculate centroids for each cluster
// // DBSCAN.prototype.getCentroids = function(dataset, labels) {
// //   var centroids = {};
// //   var clusterSums = {};
// //   var clusterCounts = {};
// //   var i, j, label;

// //   // Initialize sums and counts for each cluster
// //   for (i = 0; i < labels.length; i++) {
// //     label = labels[i];
    
// //     // Skip noise points (label -2)
// //     if (label <= 0) continue;
    
// //     if (!clusterSums[label]) {
// //       clusterSums[label] = [];
// //       clusterCounts[label] = 0;
      
// //       // Initialize sum array with zeros for each dimension
// //       for (j = 0; j < dataset[i].length; j++) {
// //         clusterSums[label][j] = 0;
// //       }
// //     }
// //   }

// //   // Sum up all points in each cluster
// //   for (i = 0; i < labels.length; i++) {
// //     label = labels[i];
    
// //     // Skip noise points
// //     if (label <= 0) continue;
    
// //     for (j = 0; j < dataset[i].length; j++) {
// //       clusterSums[label][j] += dataset[i][j];
// //     }
// //     clusterCounts[label]++;
// //   }

// //   // Calculate centroids by dividing sums by counts
// //   for (label in clusterSums) {
// //     if (clusterSums.hasOwnProperty(label)) {
// //       centroids[label] = [];
// //       for (j = 0; j < clusterSums[label].length; j++) {
// //         centroids[label][j] = clusterSums[label][j] / clusterCounts[label];
// //       }
// //     }
// //   }

// //   return centroids;
// // };

// // // Convenience method to get clusters as arrays of points
// // DBSCAN.prototype.fitPredict = function(dataset) {
// //   var labels = this.fit(dataset);
// //   var clusters = {};
// //   var noise = [];
// //   var i;

// //   for (i = 0; i < labels.length; i++) {
// //     var label = labels[i];
    
// //     if (label === -2) {
// //       noise.push({ point: dataset[i], index: i });
// //     } else {
// //       if (!clusters[label]) {
// //         clusters[label] = [];
// //       }
// //       clusters[label].push({ point: dataset[i], index: i });
// //     }
// //   }

// //   return { clusters: clusters, noise: noise, labels: labels };
// // };

// // // Enhanced method that includes centroids in the result
// // DBSCAN.prototype.fitPredictWithCentroids = function(dataset) {
// //   var labels = this.fit(dataset);
// //   var clusters = {};
// //   var noise = [];
// //   var centroids = this.getCentroids(dataset, labels);
// //   var i;

// //   for (i = 0; i < labels.length; i++) {
// //     var label = labels[i];
    
// //     if (label === -2) {
// //       noise.push({ point: dataset[i], index: i });
// //     } else {
// //       if (!clusters[label]) {
// //         clusters[label] = [];
// //       }
// //       clusters[label].push({ point: dataset[i], index: i });
// //     }
// //   }

// //   return { 
// //     clusters: clusters, 
// //     noise: noise, 
// //     labels: labels, 
// //     centroids: centroids 
// //   };
// // };

// // // Example usage:
// // function dbscan(data) {
// //   // if no data set dat to data2
// //   if (data === undefined) {
// //     data = [
// //       [8.0827, 129.55],
// //       [8.0723, 129.55],
// //       [7.5807, 132.78],
// //       [7.5693, 132.78],
// //       [7.5594, 132.78],
// //       [7.4670, 128.38],
// //       [7.4576, 128.38],
// //       [7.4488, 128.38]
// //     ];
// //   }  
// //   var dbscan1 = new DBSCAN(1.0, 1); // eps=1.0, minPts=2

// //   // // Sample 2D dataset
// //   // var data = [
// //   //   [0, 0], [1, 1], [2, 0], [8, 7], [8, 8], [25, 80],
// //   //   [0.5, 0.5], [1.5, 1.5], [9, 9], [10, 8], [7, 9]
// //   // ];

// //   // // HSQC data
// //   // var data2 = [
// //   //   [8.0827, 129.55],
// //   //   [8.0723, 129.55],
// //   //   [7.5807, 132.78],
// //   //   [7.5693, 132.78],
// //   //   [7.5594, 132.78],
// //   //   [7.4670, 128.38],
// //   //   [7.4576, 128.38],
// //   //   [7.4488, 128.38]
// //   // ];
  
// //   // Using the new method with centroids
// //   var result = dbscan1.fitPredictWithCentroids(data);
  
// //   print('Clusters:', result.clusters);
// //   print('Noise points:', result.noise);
// //   print('Labels:', result.labels);
// //   print('Centroids:', result.centroids);
// //   print('Number of clusters:', Object.keys(result.clusters).length);
// //   print('Number of centroids:', Object.keys(result.centroids).length);
  
// //   // // Or use the standalone centroid calculation
// //   // var labels = dbscan1.fit(data);
// //   // var centroids = dbscan1.getCentroids(data, labels);
// //   // print('Standalone centroids:', centroids);

// //   return result;
// // }



// // DBSCAN Clustering Implementation - ES5 Compatible
// function DBSCAN(eps, minPts) {
//   this.eps = eps || 0.5;        // Maximum distance between two points to be neighbors
//   this.minPts = minPts || 5;    // Minimum number of points to form a cluster
// }

// // Calculate Euclidean distance between two points
// DBSCAN.prototype.distance = function(point1, point2) {
//   var sum = 0;
//   for (var i = 0; i < point1.length; i++) {
//     sum += Math.pow(point1[i] - point2[i], 2);
//   }
//   return Math.sqrt(sum);
// };

// // Find all neighbors within eps distance
// DBSCAN.prototype.regionQuery = function(dataset, pointIndex) {
//   var neighbors = [];
//   var point = dataset[pointIndex];
  
//   for (var i = 0; i < dataset.length; i++) {
//     if (this.distance(point, dataset[i]) <= this.eps) {
//       neighbors.push(i);
//     }
//   }
//   return neighbors;
// };

// // Main DBSCAN clustering algorithm
// DBSCAN.prototype.fit = function(dataset) {
//   var n = dataset.length;
//   var labels = [];
//   var clusterId = 0;
//   var i, j;

//   // Initialize labels array with -1 (unclassified)
//   for (i = 0; i < n; i++) {
//     labels[i] = -1;
//   }

//   for (i = 0; i < n; i++) {
//     // Skip if point is already processed
//     if (labels[i] !== -1) continue;

//     // Find neighbors
//     var neighbors = this.regionQuery(dataset, i);

//     // Mark as noise if not enough neighbors
//     if (neighbors.length < this.minPts) {
//       labels[i] = -2; // Noise point
//       continue;
//     }

//     // Start new cluster
//     clusterId++;
//     labels[i] = clusterId;

//     // Process all neighbors
//     var seedSet = neighbors.slice(); // Copy array
//     j = 0;

//     while (j < seedSet.length) {
//       var currentPoint = seedSet[j];

//       // Change noise to border point
//       if (labels[currentPoint] === -2) {
//         labels[currentPoint] = clusterId;
//       }

//       // Skip if already processed
//       if (labels[currentPoint] !== -1) {
//         j++;
//         continue;
//       }

//       // Mark as part of cluster
//       labels[currentPoint] = clusterId;

//       // Find neighbors of current point
//       var currentNeighbors = this.regionQuery(dataset, currentPoint);

//       // If enough neighbors, add them to seed set
//       if (currentNeighbors.length >= this.minPts) {
//         for (var k = 0; k < currentNeighbors.length; k++) {
//           var neighbor = currentNeighbors[k];
//           var found = false;
          
//           // Check if neighbor already in seedSet
//           for (var m = 0; m < seedSet.length; m++) {
//             if (seedSet[m] === neighbor) {
//               found = true;
//               break;
//             }
//           }
          
//           if (!found) {
//             seedSet.push(neighbor);
//           }
//         }
//       }

//       j++;
//     }
//   }

//   return labels;
// };

// // Calculate centroids for each cluster
// DBSCAN.prototype.getCentroids = function(dataset, labels) {
//   var centroids = {};
//   var clusterSums = {};
//   var clusterCounts = {};
//   var i, j, label;

//   // Initialize sums and counts for each cluster
//   for (i = 0; i < labels.length; i++) {
//     label = labels[i];
    
//     // Skip noise points (label -2)
//     if (label <= 0) continue;
    
//     if (!clusterSums[label]) {
//       clusterSums[label] = [];
//       clusterCounts[label] = 0;
      
//       // Initialize sum array with zeros for each dimension
//       for (j = 0; j < dataset[i].length; j++) {
//         clusterSums[label][j] = 0;
//       }
//     }
//   }

//   // Sum up all points in each cluster
//   for (i = 0; i < labels.length; i++) {
//     label = labels[i];
    
//     // Skip noise points
//     if (label <= 0) continue;
    
//     for (j = 0; j < dataset[i].length; j++) {
//       clusterSums[label][j] += dataset[i][j];
//     }
//     clusterCounts[label]++;
//   }

//   // Calculate centroids by dividing sums by counts
//   for (label in clusterSums) {
//     if (clusterSums.hasOwnProperty(label)) {
//       centroids[label] = [];
//       for (j = 0; j < clusterSums[label].length; j++) {
//         centroids[label][j] = clusterSums[label][j] / clusterCounts[label];
//       }
//     }
//   }

//   return centroids;
// };

// // Convenience method to get clusters as arrays of points
// DBSCAN.prototype.fitPredict = function(dataset) {
//   var labels = this.fit(dataset);
//   var clusters = {};
//   var noise = [];
//   var i;

//   for (i = 0; i < labels.length; i++) {
//     var label = labels[i];
    
//     if (label === -2) {
//       noise.push({ point: dataset[i], index: i });
//     } else {
//       if (!clusters[label]) {
//         clusters[label] = [];
//       }
//       clusters[label].push({ point: dataset[i], index: i });
//     }
//   }

//   return { clusters: clusters, noise: noise, labels: labels };
// };

// // Standard scaler function
// DBSCAN.prototype.standardScale = function(data) {
//     if (data.length === 0) return { scaledData: data, scaleParams: null };
    
//     var scaled = [];
//     var means = [];
//     var stds = [];
//     var i, dim;
    
//     // Calculate means for each dimension
//     for (dim = 0; dim < data[0].length; dim++) {
//         var sum = 0;
//         for (i = 0; i < data.length; i++) {
//             sum += data[i][dim];
//         }
//         means[dim] = sum / data.length;
//     }
    
//     // Calculate standard deviations
//     for (dim = 0; dim < data[0].length; dim++) {
//         var sumSquaredDiffs = 0;
//         for (i = 0; i < data.length; i++) {
//             sumSquaredDiffs += Math.pow(data[i][dim] - means[dim], 2);
//         }
//         stds[dim] = Math.sqrt(sumSquaredDiffs / data.length);
//         // Avoid division by zero
//         if (stds[dim] === 0) stds[dim] = 1;
//     }
    
//     // Scale the data
//     for (i = 0; i < data.length; i++) {
//         var scaledPoint = [];
//         for (dim = 0; dim < data[i].length; dim++) {
//             scaledPoint[dim] = (data[i][dim] - means[dim]) / stds[dim];
//         }
//         scaled.push(scaledPoint);
//     }
    
//     return {
//         scaledData: scaled,
//         scaleParams: { means: means, stds: stds }
//     };
// };

// // Inverse transform to get centroids back to original scale
// DBSCAN.prototype.inverseTransform = function(centroids, scaleParams) {
//     if (!scaleParams) return centroids;
    
//     var originalCentroids = {};
//     var clusterId, dim;
    
//     for (clusterId in centroids) {
//         if (centroids.hasOwnProperty(clusterId)) {
//             originalCentroids[clusterId] = [];
//             for (dim = 0; dim < centroids[clusterId].length; dim++) {
//                 originalCentroids[clusterId][dim] = 
//                     (centroids[clusterId][dim] * scaleParams.stds[dim]) + scaleParams.means[dim];
//             }
//         }
//     }
    
//     return originalCentroids;
// };

// // Enhanced method that includes centroids in the result
// DBSCAN.prototype.fitPredictWithCentroids = function(dataset, autoScale) {
//     autoScale = autoScale !== undefined ? autoScale : true; // Default to true
    
//     var originalData = dataset;
//     var scaleParams = null;
//     var processedData = dataset;
    
//     // Scale data if requested
//     if (autoScale) {
//         var scaleResult = this.standardScale(dataset);
//         processedData = scaleResult.scaledData;
//         scaleParams = scaleResult.scaleParams;
//     }
    
//     var labels = this.fit(processedData);
//     var clusters = {};
//     var noise = [];
//     var scaledCentroids = this.getCentroids(processedData, labels);
//     var centroids = this.inverseTransform(scaledCentroids, scaleParams);
//     var i;

//     for (i = 0; i < labels.length; i++) {
//         var label = labels[i];
        
//         if (label === -2) {
//             noise.push({ point: originalData[i], index: i });
//         } else {
//             if (!clusters[label]) {
//                 clusters[label] = [];
//             }
//             clusters[label].push({ point: originalData[i], index: i });
//         }
//     }

//     return { 
//         clusters: clusters, 
//         noise: noise, 
//         labels: labels, 
//         centroids: centroids,
//         scaledData: processedData,
//         scaleParams: scaleParams
//     };
// };

// // Example usage:
// function dbscan(data, eps, minPts, autoScale) {

//   var epsValue = eps || 1.0; // Default eps value if not provided
//   var minPtsValue = minPts || 1; // Default minPts value if not provided
//   var dbscan1 = new DBSCAN(epsValue, minPtsValue); // eps=1.0, minPts=1
//   var autoScaleValue = autoScale !== undefined ? autoScale : true; // Default to true

//   // Default HSQC data if no data provided
//   var defaultData = [
//     [8.0827, 129.55],
//     [8.0723, 129.55],
//     [7.5807, 132.78],
//     [7.5693, 132.78],
//     [7.5594, 132.78],
//     [7.4670, 128.38],
//     [7.4576, 128.38],
//     [7.4488, 128.38]
//   ];
  
//   // Use provided data or default to HSQC data
//   var inputData = data || defaultData;
  
//   print('=== DBSCAN with Auto-Scaling (Recommended) ===');
//   // Using auto-scaling (default behavior)
//   var result1 = dbscan1.fitPredictWithCentroids(inputData, autoScaleValue);
  
//   print('Clusters:', result1.clusters);
//   print('Noise points:', result1.noise);
//   print('Labels:', result1.labels);
//   print('Centroids (original scale):', result1.centroids);
//   print('Number of clusters:', Object.keys(result1.clusters).length);
  
//   // print('\n=== DBSCAN without Auto-Scaling ===');
//   // // Using without scaling
//   // var result2 = dbscan1.fitPredictWithCentroids(inputData, false);
  
//   // print('Clusters:', result2.clusters);
//   // print('Centroids (no scaling):', result2.centroids);
//   // print('Number of clusters:', Object.keys(result2.clusters).length);
  
//   // print('\n=== Scale Parameters Used ===');
//   // if (result1.scaleParams) {
//   //   print('Means:', result1.scaleParams.means);
//   //   print('Standard deviations:', result1.scaleParams.stds);
//   // }
//   return result1;
// }

// DBSCAN Clustering Implementation - ES5 Compatible
function DBSCAN(eps, minPts) {
  this.eps = eps || 1;        // Maximum distance between two points to be neighbors
  this.minPts = minPts || 1;    // Minimum number of points to form a cluster
}

// Calculate Euclidean distance between two points
DBSCAN.prototype.distance = function(point1, point2) {
  var sum = 0;
  for (var i = 0; i < point1.length; i++) {
    sum += Math.pow(point1[i] - point2[i], 2);
  }
  return Math.sqrt(sum);
};

// Find all neighbors within eps distance
DBSCAN.prototype.regionQuery = function(dataset, pointIndex) {
  var neighbors = [];
  var point = dataset[pointIndex];
  
  for (var i = 0; i < dataset.length; i++) {
    if (this.distance(point, dataset[i]) <= this.eps) {
      neighbors.push(i);
    }
  }
  return neighbors;
};

// Main DBSCAN clustering algorithm
DBSCAN.prototype.fit = function(dataset) {
  var n = dataset.length;
  var labels = [];
  var clusterId = 0;
  var i, j;

  // Initialize labels array with -1 (unclassified)
  for (i = 0; i < n; i++) {
    labels[i] = -1;
  }

  for (i = 0; i < n; i++) {
    // Skip if point is already processed
    if (labels[i] !== -1) continue;

    // Find neighbors
    var neighbors = this.regionQuery(dataset, i);

    // Mark as noise if not enough neighbors
    if (neighbors.length < this.minPts) {
      labels[i] = -2; // Noise point
      continue;
    }

    // Start new cluster
    clusterId++;
    labels[i] = clusterId;

    // Process all neighbors
    var seedSet = neighbors.slice(); // Copy array
    j = 0;

    while (j < seedSet.length) {
      var currentPoint = seedSet[j];

      // Change noise to border point
      if (labels[currentPoint] === -2) {
        labels[currentPoint] = clusterId;
      }

      // Skip if already processed
      if (labels[currentPoint] !== -1) {
        j++;
        continue;
      }

      // Mark as part of cluster
      labels[currentPoint] = clusterId;

      // Find neighbors of current point
      var currentNeighbors = this.regionQuery(dataset, currentPoint);

      // If enough neighbors, add them to seed set
      if (currentNeighbors.length >= this.minPts) {
        for (var k = 0; k < currentNeighbors.length; k++) {
          var neighbor = currentNeighbors[k];
          var found = false;
          
          // Check if neighbor already in seedSet
          for (var m = 0; m < seedSet.length; m++) {
            if (seedSet[m] === neighbor) {
              found = true;
              break;
            }
          }
          
          if (!found) {
            seedSet.push(neighbor);
          }
        }
      }

      j++;
    }
  }

  return labels;
};

// Calculate centroids for each cluster
DBSCAN.prototype.getCentroids = function(dataset, labels) {
  var centroids = {};
  var clusterSums = {};
  var clusterCounts = {};
  var i, j, label;

  // Initialize sums and counts for each cluster
  for (i = 0; i < labels.length; i++) {
    label = labels[i];
    
    // Skip noise points (label -2)
    if (label <= 0) continue;
    
    if (!clusterSums[label]) {
      clusterSums[label] = [];
      clusterCounts[label] = 0;
      
      // Initialize sum array with zeros for each dimension
      for (j = 0; j < dataset[i].length; j++) {
        clusterSums[label][j] = 0;
      }
    }
  }

  // Sum up all points in each cluster
  for (i = 0; i < labels.length; i++) {
    label = labels[i];
    
    // Skip noise points
    if (label <= 0) continue;
    
    for (j = 0; j < dataset[i].length; j++) {
      clusterSums[label][j] += dataset[i][j];
    }
    clusterCounts[label]++;
  }

  // Calculate centroids by dividing sums by counts
  for (label in clusterSums) {
    if (clusterSums.hasOwnProperty(label)) {
      centroids[label] = [];
      for (j = 0; j < clusterSums[label].length; j++) {
        centroids[label][j] = clusterSums[label][j] / clusterCounts[label];
      }
    }
  }

  return centroids;
};

// Convenience method to get clusters as arrays of points
DBSCAN.prototype.fitPredict = function(dataset) {
  var labels = this.fit(dataset);
  var clusters = {};
  var noise = [];
  var i;

  for (i = 0; i < labels.length; i++) {
    var label = labels[i];
    
    if (label === -2) {
      noise.push({ point: dataset[i], index: i });
    } else {
      if (!clusters[label]) {
        clusters[label] = [];
      }
      clusters[label].push({ point: dataset[i], index: i });
    }
  }

  return { clusters: clusters, noise: noise, labels: labels };
};

// Standard scaler function
DBSCAN.prototype.standardScale = function(data) {
    if (data.length === 0) return { scaledData: data, scaleParams: null };
    
    var scaled = [];
    var means = [];
    var stds = [];
    var i, dim;
    
    // Calculate means for each dimension
    for (dim = 0; dim < data[0].length; dim++) {
        var sum = 0;
        for (i = 0; i < data.length; i++) {
            sum += data[i][dim];
        }
        means[dim] = sum / data.length;
    }
    
    // Calculate standard deviations
    for (dim = 0; dim < data[0].length; dim++) {
        var sumSquaredDiffs = 0;
        for (i = 0; i < data.length; i++) {
            sumSquaredDiffs += Math.pow(data[i][dim] - means[dim], 2);
        }
        stds[dim] = Math.sqrt(sumSquaredDiffs / data.length);
        // Avoid division by zero
        if (stds[dim] === 0) stds[dim] = 1;
    }
    
    // Scale the data
    for (i = 0; i < data.length; i++) {
        var scaledPoint = [];
        for (dim = 0; dim < data[i].length; dim++) {
            scaledPoint[dim] = (data[i][dim] - means[dim]) / stds[dim];
        }
        scaled.push(scaledPoint);
    }
    
    return {
        scaledData: scaled,
        scaleParams: { means: means, stds: stds, type: 'statistical' }
    };
};

// NMR-specific scaling based on expected chemical shift ranges
DBSCAN.prototype.nmrScale = function(data, protonRange, carbonRange) {
    if (data.length === 0) return { scaledData: data, scaleParams: null };
    
    // Default NMR ranges if not provided
    protonRange = protonRange || 14;   // 1H typical range: 0-14 ppm
    carbonRange = carbonRange || 200;  // 13C typical range: 0-200 ppm
    
    var scaled = [];
    var ranges = [protonRange, carbonRange];
    var i, dim;
    
    // Scale each point
    for (i = 0; i < data.length; i++) {
        var scaledPoint = [];
        for (dim = 0; dim < data[i].length && dim < ranges.length; dim++) {
            // Normalize to 0-1 range based on expected chemical shift range
            scaledPoint[dim] = data[i][dim] / ranges[dim];
        }
        // Handle additional dimensions beyond proton/carbon if they exist
        for (dim = ranges.length; dim < data[i].length; dim++) {
            scaledPoint[dim] = data[i][dim]; // Keep as-is for unknown dimensions
        }
        scaled.push(scaledPoint);
    }
    
    return {
        scaledData: scaled,
        scaleParams: { 
            ranges: ranges, 
            protonRange: protonRange, 
            carbonRange: carbonRange,
            type: 'nmr' 
        }
    };
};

// Inverse transform to get centroids back to original scale
DBSCAN.prototype.inverseTransform = function(centroids, scaleParams) {
    if (!scaleParams) return centroids;
    
    var originalCentroids = {};
    var clusterId, dim;
    
    if (scaleParams.type === 'statistical') {
        // Statistical scaling inverse transform
        for (clusterId in centroids) {
            if (centroids.hasOwnProperty(clusterId)) {
                originalCentroids[clusterId] = [];
                for (dim = 0; dim < centroids[clusterId].length; dim++) {
                    originalCentroids[clusterId][dim] = 
                        (centroids[clusterId][dim] * scaleParams.stds[dim]) + scaleParams.means[dim];
                }
            }
        }
    } else if (scaleParams.type === 'nmr') {
        // NMR scaling inverse transform
        for (clusterId in centroids) {
            if (centroids.hasOwnProperty(clusterId)) {
                originalCentroids[clusterId] = [];
                for (dim = 0; dim < centroids[clusterId].length; dim++) {
                    if (dim < scaleParams.ranges.length) {
                        // Scale back by multiplying with the range
                        originalCentroids[clusterId][dim] = 
                            centroids[clusterId][dim] * scaleParams.ranges[dim];
                    } else {
                        // Keep as-is for dimensions beyond proton/carbon
                        originalCentroids[clusterId][dim] = centroids[clusterId][dim];
                    }
                }
            }
        }
    }
    
    return originalCentroids;
};

// Enhanced method that includes centroids in the result
DBSCAN.prototype.fitPredictWithCentroids = function(dataset, scalingMethod, protonRange, carbonRange) {
    // scalingMethod can be: 'auto' (default), 'nmr', 'statistical', or false/null (no scaling)
    scalingMethod = scalingMethod !== undefined ? scalingMethod : 'auto';
    
    var originalData = dataset;
    var scaleParams = null;
    var processedData = dataset;
    
    // Determine scaling method
    if (scalingMethod === 'auto') {
        // Use NMR scaling for small datasets (< 20 points), statistical for larger
        scalingMethod = dataset.length < 20 ? 'nmr' : 'statistical';
    }
    
    // Scale data based on method
    if (scalingMethod === 'nmr') {
        var scaleResult = this.nmrScale(dataset, protonRange, carbonRange);
        processedData = scaleResult.scaledData;
        scaleParams = scaleResult.scaleParams;
    } else if (scalingMethod === 'statistical') {
        var scaleResult = this.standardScale(dataset);
        processedData = scaleResult.scaledData;
        scaleParams = scaleResult.scaleParams;
    }
    // If scalingMethod is false/null or unrecognized, use original data
    
    var labels = this.fit(processedData);
    var clusters = {};
    var noise = [];
    var scaledCentroids = this.getCentroids(processedData, labels);
    var centroids = this.inverseTransform(scaledCentroids, scaleParams);
    var i;

    for (i = 0; i < labels.length; i++) {
        var label = labels[i];
        
        if (label === -2) {
            noise.push({ point: originalData[i], index: i });
        } else {
            if (!clusters[label]) {
                clusters[label] = [];
            }
            clusters[label].push({ point: originalData[i], index: i });
        }
    }

    return { 
        clusters: clusters, 
        noise: noise, 
        labels: labels, 
        centroids: centroids,
        scaledData: processedData,
        scaleParams: scaleParams,
        scalingMethod: scalingMethod
    };
};



// Example usage:
function dbscan(data, eps, minPts) {
  var eps = eps || 3.0; // Default eps value if not provided
  var minPts = minPts || 1; // Default minPts value if not provided
  var dbscan1 = new DBSCAN(eps, minPts);

  // Default HSQC data if no data provided
  var defaultData = [
    [8.0827, 129.55],
    [8.0723, 129.55],
    [7.5807, 132.78],
    [7.5693, 132.78],
    [7.5594, 132.78],
    [7.4670, 128.38],
    [7.4576, 128.38],
    [7.4488, 128.38]
  ];
  
  // Use provided data or default to HSQC data
  var inputData = data || defaultData;
  
  print('=== DBSCAN with Auto-Scaling (NMR-aware for small datasets) ===');
  // Auto-scaling: uses NMR scaling for < 20 points, statistical for larger datasets
  var result1 = dbscan1.fitPredictWithCentroids(inputData, 'auto');
  
  print('Scaling method used:', result1.scalingMethod);
  print('Clusters:', result1.clusters);
  print('Noise points:', result1.noise);
  print('Labels:', result1.labels);
  print('Centroids (original scale):', result1.centroids);
  print('Number of clusters:', Object.keys(result1.clusters).length);

  return result1;
  
}