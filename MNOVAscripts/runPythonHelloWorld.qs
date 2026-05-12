function runPythonHelloWorld() {
    var proc = new Process();
    proc.start('/Users/vsmw51/opt/anaconda3/envs/py313/bin/python3', ['-c', 'print("Hello, World!")']);
    
    if (!proc.waitForFinished(5000)) {
        print('Process timed out or failed to finish');
        return null;
    }
    
    var output = proc.allStdOutput;
    var errors = proc.allErrorOutput;
    
    if (errors) {
        print('Errors: ' + errors);
    }
    
    print(output);
    return output;
}

runPythonHelloWorld();