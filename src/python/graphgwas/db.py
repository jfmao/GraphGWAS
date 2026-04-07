"""Neo4j connection wrapper. Follows GraphMana's GraphManaConnection pattern."""

from neo4j import GraphDatabase


class _EagerResult:
    """Lightweight wrapper holding materialized records from a closed session."""

    def __init__(self, records: list, summary=None):
        self._records = records
        self._summary = summary

    def single(self):
        if len(self._records) == 1:
            return self._records[0]
        if len(self._records) == 0:
            return None
        raise ValueError(f"Expected exactly one record, got {len(self._records)}")

    def __iter__(self):
        return iter(self._records)

    def __len__(self):
        return len(self._records)

    def data(self):
        return [dict(r) for r in self._records]


class GraphGWASConnection:
    """Thin wrapper around neo4j.GraphDatabase.driver().

    Usage::

        with GraphGWASConnection(uri, user, password) as conn:
            result = conn.execute_read("MATCH (n) RETURN count(n) AS c")
    """

    def __init__(self, uri: str, user: str, password: str, database: str | None = None):
        self._uri = uri
        self._user = user
        self._password = password
        self._database = database
        self._driver = None

    def __enter__(self):
        self._driver = GraphDatabase.driver(
            self._uri, auth=(self._user, self._password),
            max_connection_pool_size=50,
            connection_acquisition_timeout=60,
        )
        self._driver.verify_connectivity()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._driver:
            self._driver.close()
        return False

    def execute_read(self, query: str, parameters: dict | None = None):
        with self._driver.session(database=self._database) as session:
            result = session.run(query, parameters or {})
            return _EagerResult(list(result), result.consume())

    def execute_write(self, query: str, parameters: dict | None = None):
        with self._driver.session(database=self._database) as session:
            result = session.run(query, parameters or {})
            return _EagerResult(list(result), result.consume())

    def execute_write_tx(self, tx_func, **kwargs):
        with self._driver.session(database=self._database) as session:
            return session.execute_write(tx_func, **kwargs)

    def execute_read_tx(self, tx_func, **kwargs):
        with self._driver.session(database=self._database) as session:
            return session.execute_read(tx_func, **kwargs)

    @property
    def driver(self):
        return self._driver
