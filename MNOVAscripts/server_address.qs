// <GUI menuname="server_address"  />

function server_address() {
    // return "http://localhost:8000/";
    // On macOS, "localhost" can resolve to the IPv6 loopback (::1) ahead
    // of 127.0.0.1, and Flask's dev server here only binds to IPv4 --
    // giving an immediate "Connection refused" even though the server is
    // running. Using the literal IPv4 address sidesteps that resolution
    // order entirely.
    // return "http://127.0.0.1:5000/";
    // return "http://simplenmr.pythonanywhere.com/";
    return "http://test-simplenmr.pythonanywhere.com/";
    // return"http://simplenmr.awh.durham.ac.uk/";
}
