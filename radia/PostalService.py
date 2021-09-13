from policy.LetterCA import LetterCA
import logging
import sys
import networkx as nx

class PostalService:
    def __init__(self):
        self.CAGraph = nx.DiGraph()

        self.useful_objects = 0
        self.size_of_transmit = 0

        self.useful_objects_hist = [] # history values for plotting
        self.size_of_transmit_hist = [] # history values for plotting
        self.useful_percentage_hist = []

    def resetCAGraph(self):
        self.CAGraph.clear() 

    def postCA(self, letter):
        # print('> PS: Letter from ' + letter.origin + ' for ' + str(len(letter.recipient)) + ' recipients >')
        for recipient in letter.recipient: 
        #     glob = [obj for obj in self.globs if obj.id==recipient]
        #     if len(glob) == 0:
        #         sys.exit('Error: a recipient is not available in the global object list')
        #     elif len(glob) > 1:
        #         sys.exit('There is more than two objects in the global object list with the same ID')
            # ob = glob[0] # the true global objects
            # ob.acquaintances.append(letter.origin) if letter.origin not in ob.acquaintances else ob.acquaintances # actually this is not necessary since acquaintances have been obtained and position can be directly obtained from global object list
            self.CAGraph.add_edges_from([(letter.origin, recipient)])
    
    def resetCPFloodGraph(self): # should be called before CP exchange round
        self.useful_objects = 0
        self.size_of_transmit = 0

    def postCPFloodUseful(self, useful_objects):
        self.useful_objects += useful_objects

    def postCPFloodTransmitSize(self, size_of_transmit):
        self.size_of_transmit += size_of_transmit

    def postCPFloodSaveHistory(self):
        self.useful_objects_hist.append(self.useful_objects)
        self.size_of_transmit_hist.append(self.size_of_transmit)
        if self.size_of_transmit == 0:
            useful_percentage = 1
        else:
            useful_percentage = self.useful_objects/self.size_of_transmit
        self.useful_percentage_hist.append(useful_percentage)