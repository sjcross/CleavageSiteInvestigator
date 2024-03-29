class ErrorStore():
    def __init__(self):
        self._store = {'NO_CASS':0,'CASS_MISMATCH':0,'NO_CONSENSUS':0,'CONSENSUS_MISMATCH':0,'MAX_GAP':0,'NO_MIDPOINT':0}

    def get_store(self):
        return self._store

    def cassette_not_found_in_consensus(self):
        self._increment_counter('NO_CASS')

    def cassette_ends_mismatch(self):
        self._increment_counter('CASS_MISMATCH')

    def consensus_not_found_in_reference(self):
        self._increment_counter('NO_CONSENSUS')

    def consensus_ends_mismatch(self):
        self._increment_counter('CONSENSUS_MISMATCH')

    def max_gap_exceeded(self):
        self._increment_counter('MAX_GAP')

    def midpoint_not_found(self):
        self._increment_counter('NO_MIDPOINT')

    def print_counts(self, offset=""):
        print(f"{offset}Cassette not found in consensus:  {self._store['NO_CASS']}")
        print(f"{offset}Cassette ends RC mismatch:   {self._store['CASS_MISMATCH']}")
        print(f"{offset}Consensus not found in reference: {self._store['NO_CONSENSUS']}")
        print(f"{offset}Consensus ends RC mismatch:       {self._store['CONSENSUS_MISMATCH']}")
        print(f"{offset}Maximum gap exceeded:        {self._store['MAX_GAP']}")
        print(f"{offset}Midpoint not found:          {self._store['NO_MIDPOINT']}")
        print("\n")

    def _increment_counter(self, label):
        self._store[label] = self._store[label] + 1
