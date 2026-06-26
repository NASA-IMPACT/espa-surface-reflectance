import os
import requests
import logging


# Determine the logging level. Default is WARNING.
espa_log_level = os.environ.get('ESPA_LOG_LEVEL')
if espa_log_level == 'DEBUG':
    log_level = logging.DEBUG
elif espa_log_level == 'INFO':
    log_level = logging.INFO
else:
    log_level = logging.WARNING
logging.getLogger('requests').setLevel(log_level)


class APIException(Exception):
    """
    Handle exceptions thrown by the APIServer class
    """
    pass


class APIServer(object):
    """
    Provide a more straightforward way of handling API calls
    without changing the cron jobs significantly
    """
    def __init__(self, base_url):
        self.base = base_url

    def request(self, method, resource=None, status=None, **kwargs):
        """
        Make a call into the API
        Args:
            method: HTTP method to use
            resource: API resource to touch
        Returns: response and status code
        """
        valid_methods = ('get', 'put', 'delete', 'head', 'options', 'post')

        if method not in valid_methods:
            raise APIException('Invalid method {}'.format(method))

        if resource and resource[0] == '/':
            url = '{}{}'.format(self.base, resource)
        elif resource:
            url = '{}/{}'.format(self.base, resource)
        else:
            url = self.base

        try:
            resp = requests.request(method, url, **kwargs)
        except requests.RequestException as e:
            raise APIException(e)

        if status and resp.status_code != status:
            self._unexpected_status(resp.status_code, url)

        return resp.json(), resp.status_code

    def get_configuration(self, key):
        """
        Retrieve a configuration value
        Args:
            key: configuration key
        Returns: value if it exists, otherwise None
        """
        config_url = '/configuration/{}'.format(key)

        resp, status = self.request('get', config_url, status=200)

        if key in list(resp.keys()):
            return resp[key]

    @staticmethod
    def _unexpected_status(code, url):
        """
        Throw exception for an unhandled http status
        Args:
            code: http status that was received
            url: URL that was used
        """
        raise Exception('Received unexpected status code: {}\n'
                        'for URL: {}'.format(code, url))

    def test_connection(self):
        """
        Tests the base URL for the class
        Returns: True if 200 status received, else False
        """
        resp, status = self.request('get')

        if status == 200:
            return True

        return False


def api_connect(url):
    """
    Simple lead in method for using the API connection class
    Args:
        url: base URL to connect to
    Returns: initialized APIServer object if successful connection
             else None
    """
    api = APIServer(url)

    if not api.test_connection():
        return None

    return api
