import React ,{Component} from 'react'
import Test from '../components/test'
import { connect } from 'react-redux'
import { bindActionCreators } from 'redux'
import * as testActions from '../actions/testAction'

function mapStateToProps(store) {
	return{
		data:store.test.data
	}
}

function mapDispatchToProps(dispatch) {
  return {
    actions: bindActionCreators(testActions, dispatch)
  }
}


let TestContainer=connect(mapStateToProps,mapDispatchToProps)(Test)

export default TestContainer