from radia import angutil as ang
import numpy as np
from radia import v2x as v2x
from radia import LetterCA
import sys

class GlobalMapObject:
    def __init__(self, id, pos, spd, acc, orientation):
        self.id = id
        self.pos = pos
        self.posx = np.array(list(pos))
        self.spd = spd
        self.acc = acc
        self.orientation = ang.roundang(orientation * 180 / np.pi)
        self.neighbors = [] # detectable by radar 
        self.acquaintances = []  # communicable via V2X
        self.informs_candidate = [] # objects received but yet to be transmitted (due to loop synchronization, make sure all vehicles are done receiving first)
        self.informs = [] # objects received and ready for transmission
        self.memory = {} # memory about past vehicle trajectories
        self.otqd_entropy = {} # likelihood of observed trajectories being nominal

    def __eq__(self, other):
        if not isinstance(other, GlobalMapObject):
            return NotImplemented
        return self.id==other.id

    def newTransmissionRound(self):
        self.neighbors.clear()
        self.acquaintances.clear()
        self.informs = self.informs_candidate
        self.informs_candidate = []

    def addNeighbor(self, id):
        self.neighbors.append(id) if id not in self.neighbors else self.neighbors

    def addAcquaintance(self, id):
        self.acquaintances.append(id) if id not in self.acquaintances else self.acquaintances

    def transmitCAM(self, postal_service):
        l = LetterCA(self.acquaintances, self.posx, self.id, self.posx)
        postal_service.postCA(l) # update the comm graph

    def getPos(self, globs, id):
        obj = [x for x in globs if x.id==id]
        if len(obj)==1:
            return obj
        else:
            sys.exit('Requested ID ' + id + ' is invalid')
    
    def transmitCPFlood(self, globs, postal_service):
        # filter out relayed distant objects 
        for obj in self.informs:
            if np.linalg.norm(obj.posx - self.posx) > v2x.DROP_RANGE:
                self.informs.remove(obj)
        to_transmit = set(self.neighbors).union(set(self.informs)) # all objects known to self
        postal_service.postCPFloodTransmitSize(len(to_transmit))
        for obj in globs: 
            if obj.id in self.acquaintances:
                client_informed = set(obj.informs_candidate) # all objects transmitted to client earlier (by other stations)
                pre_union_card = len(client_informed)
                client_informed = client_informed.union(to_transmit) # union with all objects known to self
                post_union_card = len(client_informed)
                if post_union_card - pre_union_card > 0:
                    print('Increase in ', obj.id, ' perceptive field: ', str(post_union_card-pre_union_card), ' objects were added from ', self.id)
                obj.informs_candidate = list(client_informed) 

    def update_statistics_before_CP(self):
        self.known_len = len(set(self.neighbors).union(set(self.informs)).union(set(self.informs_candidate)))

    def update_statistics_after_CP(self, postal_service):
        postal_service.postCPFloodUseful(len(set(self.neighbors).union(set(self.informs)).union(set(self.informs_candidate))) - self.known_len)