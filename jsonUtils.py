def extract_hostname(json_data:dict) -> str:
    """
    Extracts the host name from the given JSON data.

    Args:
        json_data (dict): The JSON data containing the host name.

    Returns:
        str: The extracted hostname.
    """
    try:
        # Assuming the host ID is located at json_data['hostId']
        print("json_data['hostname']:", json_data['hostname'])
        return json_data['hostname']["data"]["0"]
    except KeyError:
        return ""