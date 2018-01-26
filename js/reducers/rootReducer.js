import { combineReducers } from 'redux'
import testReducer from './test/testReducer'


const rootReducer = combineReducers({
  test: testReducer,
  //user: userReducer
})

export default rootReducer
