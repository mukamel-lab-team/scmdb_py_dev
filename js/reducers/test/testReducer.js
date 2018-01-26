
export default function testReducer(testData={},action) {
	switch(action.type){
         case 'TEST_ACTION':
            //do any mutation in store here
            return Object.assign({},testData);
            break

		default :
		  return testData;
		  break;
	}
}