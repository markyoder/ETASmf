
class test_class(object):
	def __init__(self):
		import matplotlib.pyplot as plt
		self.plt=plt
	#
	def dofig(self, fnum):
		self.f=self.plt.figure(fnum)
		
