

// DBSCAN Clustering Implementation - ES5 Compatible
function DBSCAN(eps, minPts) {
  this.eps = eps || 0.5;        // Maximum distance between two points to be neighbors
  this.minPts = minPts || 5;    // Minimum number of points to form a cluster
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

// Enhanced method that includes centroids in the result
DBSCAN.prototype.fitPredictWithCentroids = function(dataset) {
  var labels = this.fit(dataset);
  var clusters = {};
  var noise = [];
  var centroids = this.getCentroids(dataset, labels);
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

  return { 
    clusters: clusters, 
    noise: noise, 
    labels: labels, 
    centroids: centroids 
  };
};

// Example usage:
function dbscan(data) {
  // if no data set dat to data2
  if (data === undefined) {
    data = [
      [8.0827, 129.55],
      [8.0723, 129.55],
      [7.5807, 132.78],
      [7.5693, 132.78],
      [7.5594, 132.78],
      [7.4670, 128.38],
      [7.4576, 128.38],
      [7.4488, 128.38]
    ];
  }  
  var dbscan1 = new DBSCAN(1.0, 1); // eps=1.0, minPts=2

  // // Sample 2D dataset
  // var data = [
  //   [0, 0], [1, 1], [2, 0], [8, 7], [8, 8], [25, 80],
  //   [0.5, 0.5], [1.5, 1.5], [9, 9], [10, 8], [7, 9]
  // ];

  // // HSQC data
  // var data2 = [
  //   [8.0827, 129.55],
  //   [8.0723, 129.55],
  //   [7.5807, 132.78],
  //   [7.5693, 132.78],
  //   [7.5594, 132.78],
  //   [7.4670, 128.38],
  //   [7.4576, 128.38],
  //   [7.4488, 128.38]
  // ];
  
  // Using the new method with centroids
  var result = dbscan1.fitPredictWithCentroids(data);
  
  print('Clusters:', result.clusters);
  print('Noise points:', result.noise);
  print('Labels:', result.labels);
  print('Centroids:', result.centroids);
  print('Number of clusters:', Object.keys(result.clusters).length);
  print('Number of centroids:', Object.keys(result.centroids).length);
  
  // // Or use the standalone centroid calculation
  // var labels = dbscan1.fit(data);
  // var centroids = dbscan1.getCentroids(data, labels);
  // print('Standalone centroids:', centroids);

  return result;
}

