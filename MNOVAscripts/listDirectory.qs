function listDirectory(path) {
    path = path || '.';
    
    var proc = new Process();
    proc.start('ls', [path]);
    
    if (!proc.waitForFinished(5000)) {
        print('Process timed out or failed to finish');
        return null;
    }
    
    var output = proc.allStdOutput;
    var errors = proc.allErrorOutput;
    
    if (errors) {
        print('Errors: ' + errors);
    }
    
    return output;
}

// Usage
var result = listDirectory();
print(result);