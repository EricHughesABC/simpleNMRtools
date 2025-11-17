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

// Example usage:
var dbscan = new DBSCAN(1.0, 2); // eps=1.0, minPts=3

// Sample 2D dataset
var data = [
  [0, 0], [1, 1], [2, 0], [8, 7], [8, 8], [25, 80],
  [0.5, 0.5], [1.5, 1.5], [9, 9], [10, 8], [7, 9]
];

var result = dbscan.fitPredict(data);

print('Clusters:', result.clusters);
print('Noise points:', result.noise);
print('Labels:', result.labels);