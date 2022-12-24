import bisect

lower = bisect.bisect_right()


#
# db = GeneDetector('C:/Users/keren/PycharmProjects/project/Homo_sapiens.GRCh37.87.gff3')
# genes = db.test()
# genes.to_csv("genes", sep='\t', encoding='utf-8')
#
# #inside class
#
# def test(self):
#     return self.genes_db
#
# def tester():
#     print(db.inside_gene(1, 69091, 70008) == True) #exact
#     print(db.inside_gene(1, 1109265, 1133215) == True) #inside
#     print(db.inside_gene(1, 16865516, 16866070) == False) #out left
#     print(db.inside_gene(1, 19991780, 20006056) == False) #out right
#     print(db.inside_gene(10, 69091, 70008) == False) # different chrom
#     print(db.inside_gene(10, 92888, 300577) == False) #overlap
#     print(db.inside_gene(15, 23810500, 23873000) == True)
#     print(db.inside_gene(19, 44088520, 44100297) == False) #out both
#     print(db.inside_gene(36, 60000, 99000) == True) #str name
#     print(db.inside_gene(36, 23516, 25354) == False)  # different str name
#     print(db.inside_gene(23, 1010, 1256) == True) #str name
#     print(db.inside_gene(51, 30333, 32000) == True) ###############
#     print(db.inside_gene(51, 36277,	50281) == True)
#     print(db.inside_gene(26,1,1000) == False)
#     print(db.inside_gene(26,410000,420000) == False)
#     print(db.inside_gene(29,49252,88375) == True)
#     print(db.inside_gene(29,86545,86855) == True)
#     print(db.inside_gene(28,53590,56132) == True)
#     print(db.inside_gene(33,43717,74342) == True)
#     print(db.inside_gene(77,3905,22432) == True)
#     print(db.inside_gene(77,26348,27637) == True)
#     print(db.inside_gene(54,44537,47211) == True)
#     print(db.inside_gene(37,117169,119085) == True)
#     print(db.inside_gene(46,7891,96246) == True)
#     print(db.inside_gene(46,147649,156325) == True)
#     print(db.inside_gene(32,64252,65919) == True)
#     print(db.inside_gene(43,68442,70280) == True)
#     print(db.inside_gene(43,108007,139339) == True)
#     print(db.inside_gene(38,23516,25354) == True)
#     print(db.inside_gene(45,38792,97421) == True)
#     print(db.inside_gene(36,53831,99642) == True)
#     print(db.inside_gene(47,62571,62753) == True)
#     print(db.inside_gene(31,86430,89304) == True)
#     print(db.inside_gene(31,129789,134090) == True)
#     print(db.inside_gene(34,42781,68787) == True)
#     print(db.inside_gene(34,93511,119608) == True)
#     print(db.inside_gene(34,149787,180454) == True)
#     print(db.inside_gene(49,30530,32226) == True)
#     print(db.inside_gene(56,899,2487) == True)
#     print(db.inside_gene(58,12836,34543) == True)
#     ####### X chromosome
#     print(db.inside_gene(24, 19970000, 19980000) == True)
#     print(db.inside_gene(24, 19988417, 19988420) == False)
#     print(db.inside_gene(24, 20000000, 20135030) == False)
#     print(db.inside_gene(24, 21674444, 21676400) == True)
#     print(db.inside_gene(24, 25021811, 25034065) == True)
#     print(db.inside_gene(24, 192980, 220024) == False)
#     print(db.inside_gene(24, 14747, 15887) == False)
#     print(db.inside_gene(24, 1710486, 1721407) == True)
#     print(db.inside_gene(24, 2770781, 2771802) == True)
#     print(db.inside_gene(24, 2770781, 2771801) == True)
#     ##### 3 chromosome
#     print(db.inside_gene(3, 51812580, 51813009) == True)
#     print(db.inside_gene(3, 51812500, 51813548) == False)
#     print(db.inside_gene(3, 52232102, 52248343) == True)
#     print(db.inside_gene(3, 52288437, 52322036) == True)
#     print(db.inside_gene(3, 52321105, 52329272) == True)
#     print(db.inside_gene(3, 52350335, 52434507) == True)
#     print(db.inside_gene(3, 52435030, 52444365) == True)
#     print(db.inside_gene(3, 52444672, 52457658) == False)
#     print(db.inside_gene(3, 52467069, 52479101) == True)
#     print(db.inside_gene(3, 52485118, 52488086) == True)
#     ##### 19 chromosome
#     print(db.inside_gene(19, 107471, 111696) == True)
#     print(db.inside_gene(19, 281040, 291393) == True)
#     print(db.inside_gene(19, 1, 2) == False)
#     print(db.inside_gene(19, 2, 10) == False)
#     print(db.inside_gene(19, 1000, 1010) == False)
#     print(db.inside_gene(19, 2390, 3000) == False)
#     print(db.inside_gene(19, 1, 2000) == False)
#     print(db.inside_gene(19, 9999, 10000) == False)
#     print(db.inside_gene(19, 305573, 306467) == True)
#     print(db.inside_gene(19, 406000, 407000) == True)
#     ##### chromosome 18
#     print(db.inside_gene(18, 77941563, 77960875) == True)
#     print(db.inside_gene(18, 77915115, 78005429) == True)
#     print(db.inside_gene(18, 77905807, 77936315) == True)
#     print(db.inside_gene(18, 77866915, 77905406) == True)
#     print(db.inside_gene(18, 77866910, 77905406) == False)
#     print(db.inside_gene(18, 77794358, 77806397) == True)
#     print(db.inside_gene(18, 77854969, 77854970) == False)
#     print(db.inside_gene(18, 77732867, 77793949) == True)
#     print(db.inside_gene(18, 77724561, 77730822) == True)
#     print(db.inside_gene(18, 77960880, 77960881) == True)
#     ##### 7 chromosome
#     print(db.inside_gene(7, 26706681, 27034858) == True)
#     print(db.inside_gene(7, 27132612, 27135615) == True)
#     print(db.inside_gene(7, 38762563, 38971994) == True)
#     print(db.inside_gene(7, 39585245, 39605838) == True)
#     print(db.inside_gene(7, 65112081, 65113280) == True)
#     print(db.inside_gene(7, 65670186, 65885530) == True)
#     print(db.inside_gene(7, 66386212, 66423538) == True)
#     print(db.inside_gene(7, 77166592, 77269388) == True)
#     print(db.inside_gene(7, 90030000, 90040000) == True)
#     print(db.inside_gene(7, 90095738, 90839905) == True)
#     ##### chromosome 16
#     print(db.inside_gene(16, 1777777, 1781233) == True)
#     print(db.inside_gene(16, 1763600, 1764150) == True)
#     print(db.inside_gene(16, 1801560, 1802263) == True)
#     print(db.inside_gene(16, 1811234, 1811235) == True)
#     print(db.inside_gene(16, 1819000, 1819090) == True)
#     print(db.inside_gene(16, 1820287, 1820400) == True)
#     print(db.inside_gene(16, 1821891, 1823156) == True)
#     print(db.inside_gene(16, 2009500, 2009505) == False)
#     print(db.inside_gene(16, 2012053, 2014861) == True)
#     print(db.inside_gene(16, 2016824, 2018976) == True)
#
#
#
#
#
#
#
#
#
#
#
# start_time = time.time()
#
# tester()
#
# print("--- %s seconds ---" % (time.time() - start_time))
#
#
#
#
#
#



