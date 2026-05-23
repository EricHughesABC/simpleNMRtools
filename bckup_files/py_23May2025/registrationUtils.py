registered_hosts = ["host1", "host2", "host3", "Z49BR-HKPP041Y-42DHS-ANANTMPB"]


def is_registered(hostName):
    """
    Check if the host is registered in the system.
    :param hostName: The name of the host to check.
    :return: True if registered, False otherwise.
    """
    # Placeholder for actual registration check logic
    # In a real scenario, this might involve checking a database or an API
    return hostName in registered_hosts


def register_host(hostid, email):
    """
    Register a host in the system.
    :param hostid: The ID of the host to register.
    :param email: The email address associated with the registration.
    :return: True if registration is successful, False otherwise.
    """
    # Placeholder for actual registration logic
    # In a real scenario, this might involve sending a request to a server
    print(f"Registering host {hostid} with email {email}")

    # append the host ID to the registered hosts list
    registered_hosts.append(hostid)
    return True
